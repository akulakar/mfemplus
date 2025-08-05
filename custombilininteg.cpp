// This file containts custom integrators for the MFEM library..
// List of integrators is given below.
// 1. Elasticity
//      (a) Isotropic elasticity with E and nu.
//      (b) Fully anisotropic elasticty with 21 independent constants.
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// 2.

#include "custombilininteg.hpp"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{

    void IsotropicElasticityIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        mfem::real_t w, E, NU;

        MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
        DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
        Vector divshape(dim * dof);
#else
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);

#endif

        elmat.SetSize(dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        mfem::DenseMatrix C(str_comp, str_comp);  // Stiffness in Voigt form
        mfem::DenseMatrix B(str_comp, dof * dim); // Strain displacement matrix

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            NU = poisson_ratio->Eval(Trans, ip);
            E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            if (dim == 2)
            {
                C = 0.0;
                // Plane strain
                C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                C(2, 2) = E / (2 * (1 + NU));

                // Plane stress
                // C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                // C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                // C(2, 2) = (E * (1 - NU) / (2 * (1 - pow(NU, 2))));

                // 3D to 2D, but this is improper
                // C(0, 0) = C(1, 1) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                // C(0, 1) = C(1, 0) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                // C(2, 2) =  E / (2 * (1 + NU));

                // In 2D, we have 3 unique strain components.
                B = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
            }

            else if (dim == 3)
            {
                C = 0.0;
                C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));

                // In 3D, we have 6 unique strain components.
                B = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = gshape(spf, 2);
                    B(3, spf + 2 * dof) = gshape(spf, 1);
                    B(4, spf) = gshape(spf, 2);
                    B(4, spf + 2 * dof) = gshape(spf, 0);
                    B(5, spf) = gshape(spf, 1);
                    B(5, spf + dof) = gshape(spf, 0);
                }
            }
            mfem::DenseMatrix CB(str_comp, dof * dim);
            mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_intpt); // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w, elmat_intpt);
        }
    }

    void AnisotropicElasticityIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        mfem::real_t w;

        MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
        DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
        Vector divshape(dim * dof);
#else
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);

#endif

        elmat.SetSize(dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 3 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        mfem::DenseMatrix B(str_comp, dof * dim); // Strain displacement matrix
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            mfem::DenseMatrix C(str_comp, str_comp);
            stiffness->Eval(C, Trans, ip); // The stiffness matrix is evaluated at each integration point.

            if (dim == 2) // Plain strain
            {
                // In 2D, we have 3 unique strain components.
                B = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
            }

            else if (dim == 3)
            {
                // In 3D, we have 6 unique strain components.
                B = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = gshape(spf, 2);
                    B(3, spf + 2 * dof) = gshape(spf, 1);
                    B(4, spf) = gshape(spf, 2);
                    B(4, spf + 2 * dof) = gshape(spf, 0);
                    B(5, spf) = gshape(spf, 1);
                    B(5, spf + dof) = gshape(spf, 0);
                }
            }
            mfem::DenseMatrix CB(str_comp, dof * dim);
            mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_intpt); // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w, elmat_intpt);
        }
    }

    void IsotropicElasticityDeviatoricIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                                                        mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 4 : 6; // Deviatoric strain in 2D has 4 unique components using plane strain assumption, not 3!
        mfem::real_t w, MU;

        MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
        DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
        Vector divshape(dim * dof);
#else
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        divshape.SetSize(dim * dof);
#endif

        elmat.SetSize(dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        mfem::DenseMatrix B(str_comp, dof * dim), B_vol(str_comp, dof * dim); // Strain displacement matrix

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            MU = shear_mod->Eval(Trans, ip); // The shear modulus is evaluated at each integration point.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} 2 \mu B_dev^T B_dev. The B matrix as 3 rows in 2D and 6 rowd in 3D.
            // B_dev IS NOT the same as B. It only constructs the deviatoric part.
            if (dim == 2)
            {
                // In 2D, we are using the plane strain assumption. The deviatoric strain has 4 unique components, not 3!
                B = 0.0;
                B_vol = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(3, spf) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 1);
                    B(3, spf + dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 0);
                }

                for (int dimension = 0; dimension < dim; dimension++)
                {
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B_vol(0, spf + (dimension * dof)) = gshape(spf, dimension);
                        B_vol(1, spf + (dimension * dof)) = gshape(spf, dimension);
                        B_vol(2, spf + (dimension * dof)) = gshape(spf, dimension);
                    }
                }
            }

            else if (dim == 3)
            {
                // In 3D, we have 6 unique strain components.
                B = 0.0;
                B_vol = 0.0;

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                    B(3, spf + dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 2);
                    B(3, spf + 2 * dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 1);
                    B(4, spf) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 2);
                    B(4, spf + 2 * dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 0);
                    B(5, spf) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 1);
                    B(5, spf + dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 0);
                }
                for (int dimension = 0; dimension < dim; dimension++)
                {
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B_vol(0, spf + (dimension * dof)) = gshape(spf, dimension);
                        B_vol(1, spf + (dimension * dof)) = gshape(spf, dimension);
                        B_vol(2, spf + (dimension * dof)) = gshape(spf, dimension);
                    }
                }
            }
            mfem::DenseMatrix B_dev(B);
            B_dev.Add(-1.0 / 3, B_vol);
            mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
            mfem::MultAtB(B_dev, B_dev, elmat_intpt);
            elmat.Add(w * 2.0 * MU, elmat_intpt);
        }
    }
}
