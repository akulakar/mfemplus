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
        mfem::DenseMatrix CB(str_comp, dof * dim);
        mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);

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
            // \int_{\Omega} 2 \mu B_dev^T B_dev. The B matrix has 3 rows in 2D and 6 rowd in 3D.
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
                    B(3, spf) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 1); // This line's index is 3 for construction of deviatoric strain.
                    B(3, spf + dof) = (1.0 / pow(2.0, 0.5)) * gshape(spf, 0);

                    // There a factor of 1 / sqrt(2) and not 1 / 2 for the shear components because
                    // in Voigt form the strain vector has coefficients of sqrt(2) for shear strains.
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

                    // There a factor of 1 / sqrt(2) and not 1 / 2 for the shear components because
                    // in Voigt form the strain vector has coefficients of sqrt(2) for shear strains.
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

    void GLVectorMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans,
                                                       mfem::DenseMatrix &elmat)
    {
        int nd = el.GetDof();
        int spaceDim = Trans.GetSpaceDim();

        mfem::real_t norm;

        // If vdim is not set, set it to the space dimension
        vdim = (vdim == -1) ? spaceDim : vdim;

        elmat.SetSize(nd * vdim);
        shape.SetSize(nd);
        partelmat.SetSize(nd);
        if (VQ)
        {
            vec.SetSize(vdim);
        }
        else if (MQ)
        {
            mcoeff.SetSize(vdim);
        }

        // const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        mfem::IntegrationRule ir;
        if (el.GetGeomType() == mfem::Geometry::CUBE)
        {
            mfem::IntegrationRule *ir1D = new mfem::IntegrationRule;
            mfem::QuadratureFunctions1D::GaussLobatto(el.GetOrder() + 1, ir1D);
            ir = (mfem::IntegrationRule(*ir1D, *ir1D, *ir1D));
            delete ir1D;
        }
        else if (el.GetGeomType() == mfem::Geometry::SQUARE)
        {
            mfem::IntegrationRule *ir1D = new mfem::IntegrationRule;
            mfem::QuadratureFunctions1D::GaussLobatto(el.GetOrder() + 1, ir1D);
            ir = (mfem::IntegrationRule(*ir1D, *ir1D));
            delete ir1D;
        }
        else
        {
            int order = 2 * el.GetOrder() + Trans.OrderW() + Q_order;

            if (el.Space() == mfem::FunctionSpace::rQk)
            {
                ir = mfem::RefinedIntRules.Get(el.GetGeomType(), order);
            }
            else
            {
                ir = mfem::IntRules.Get(el.GetGeomType(), order);
            }
        }

        // const mfem::IntegrationRule *ir(&(el.GetNodes()));

        elmat = 0.0;
        for (int s = 0; s < ir.GetNPoints(); s++)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(s);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            norm = ip.weight * Trans.Weight();

            MultVVt(shape, partelmat);

            if (VQ)
            {
                VQ->Eval(vec, Trans, ip);
                for (int k = 0; k < vdim; k++)
                {
                    elmat.AddMatrix(norm * vec(k), partelmat, nd * k, nd * k);
                }
            }
            else if (MQ)
            {
                MQ->Eval(mcoeff, Trans, ip);
                for (int i = 0; i < vdim; i++)
                    for (int j = 0; j < vdim; j++)
                    {
                        elmat.AddMatrix(norm * mcoeff(i, j), partelmat, nd * i, nd * j);
                    }
            }
            else
            {
                if (Q)
                {
                    norm *= Q->Eval(Trans, ip);
                }
                partelmat *= norm;
                for (int k = 0; k < vdim; k++)
                {
                    elmat.AddMatrix(partelmat, nd * k, nd * k);
                }
            }
        }
    };

    //-----------------------------------------------------------------------------------------------------------------------
    GLMassIntegrator::GLMassIntegrator(const mfem::IntegrationRule *ir)
        : mfem::BilinearFormIntegrator(ir), Q(nullptr), maps(nullptr), geom(nullptr)
    {
        // static GLMassIntegrator::Kernels kernels;
    }

    GLMassIntegrator::GLMassIntegrator(mfem::Coefficient &q, const mfem::IntegrationRule *ir)
        : GLMassIntegrator(ir)
    {
        Q = &q;
    }

    void GLMassIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Trans,
                                                 mfem::DenseMatrix &elmat)
    {
        int nd = el.GetDof();
        // int dim = el.GetDim();
        mfem::real_t w;

#ifdef MFEM_THREAD_SAFE
        Vector shape;
#endif
        elmat.SetSize(nd);
        shape.SetSize(nd);

        mfem::IntegrationRule ir;
        if (el.GetGeomType() == mfem::Geometry::CUBE)
        {
            mfem::IntegrationRule *ir1D = new mfem::IntegrationRule;
            mfem::QuadratureFunctions1D::GaussLobatto(el.GetOrder() + 1, ir1D);
            ir = (mfem::IntegrationRule(*ir1D, *ir1D, *ir1D));
            delete ir1D;
        }
        else if (el.GetGeomType() == mfem::Geometry::SQUARE)
        {
            mfem::IntegrationRule *ir1D = new mfem::IntegrationRule;
            mfem::QuadratureFunctions1D::GaussLobatto(el.GetOrder() + 1, ir1D);
            ir = (mfem::IntegrationRule(*ir1D, *ir1D));
            delete ir1D;
        }
        else
        {
            const mfem::IntegrationRule ir = *GetIntegrationRule(el, Trans);
        }

        elmat = 0.0;
        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(i);
            Trans.SetIntPoint(&ip);

            el.CalcPhysShape(Trans, shape);

            w = Trans.Weight() * ip.weight;
            if (Q)
            {
                w *= Q->Eval(Trans, ip);
            }

            AddMult_a_VVt(w, shape, elmat);
        }
    }

    void IsotropicElasticityDamageIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        // uses spectral split of elastic strain energy.
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
        shape.SetSize(dof);

#endif
        elmat.SetSize(dof * dim);

        mfem::Array<int> eldofs(dof); // scalar for damage
        mfem::Vector eldofdamage(dof);
        eldofs = 0;
        eldofdamage = 0.0;

        int elnum = Trans.ElementNo;

        damage_fes->GetElementDofs(elnum, eldofs);

        for (int i = 0; i < eldofdamage.Size(); i++)
        {
            eldofdamage(i) = (*damage_gf)(eldofs[i]);
        }

        // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
        // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        mfem::DenseMatrix C(str_comp, str_comp);  // Stiffness in Voigt form
        mfem::DenseMatrix B(str_comp, dof * dim); // Strain displacement matrix
        mfem::DenseMatrix CB(str_comp, dof * dim);
        mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
        double damage_value, degradation_constant;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            // el.CalcShape(ip, shape);
            // k_eps_val not expected to change spatially.
            // multiply degradation constant to the final element matrix.

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);
            damage_value = mfem::InnerProduct(shape, eldofdamage);
            degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

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

            mfem::Mult(C, B, CB);                             // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_intpt);                // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w * degradation_constant, elmat_intpt); // multiplied by the degradation constant.
        }
    }

    void IsotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
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
        shape.SetSize(dof);

#endif
        elmat.SetSize(dof);

        mfem::Array<int> eldofs(dof * dim); // vector valued for displacement
        mfem::Vector eldofdisp(dof * dim);

        int elnum = Trans.ElementNo;

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        mfem::DenseMatrix C(str_comp, str_comp);   // Stiffness in Voigt form
        mfem::DenseMatrix B(str_comp, dof * dim);  // Strain displacement matrix
        mfem::DenseMatrix CB(str_comp, dof * dim); // Stiffness times strain displacement
        mfem::Vector CBu(str_comp);
        mfem::Vector Bu(str_comp);
        mfem::DenseMatrix elmat_intpt(dof, dof);
        double lambda1(0.0), lambda2(0.0), lambda3(0.0);
        double strain_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            w = ip.weight * Trans.Weight(); // Quadrature weights

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

            // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
            // This is equivalent to.
            mfem::Mult(C, B, CB);    // CB is 6 x (dof * dim)
            CB.Mult(eldofdisp, CBu); // CBu has dimension strain_comps. This is the stress vector.

            // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            if (dim == 2)
            {
                // In 2D lambda min is lambda1.
                lambda1 = (CBu(0) + CBu(1)) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
                lambda2 = (CBu(0) + CBu(1)) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
                lambda3 = lambda1 + 1.0;
            }
            if (dim == 3)
            {
                lambda1 = CBu(0) - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
                lambda2 = CBu(1) - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
                lambda3 = CBu(2) - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            }

            if (std::min({lambda1, lambda2, lambda3}) > 0)
            {
                B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
                strain_energy = mfem::InnerProduct(CBu, Bu);
            }
            else
                strain_energy = 1.0e-40;

            mfem::AddMult_a_VVt(w * strain_energy, shape, elmat); // multiplied by twice the strain energy.
        }
    }
}
