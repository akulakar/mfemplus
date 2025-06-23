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

    void ThreeDIsotropicElasticityIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
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

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            NU = poisson_ratio->Eval(Trans, ip);
            E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.

            mfem::DenseMatrix C; // Stiffness in Voigt notation.
            C.SetSize(6, 6);
            C = 0.0;

            C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
            C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
            C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));

            // Here we want to use voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B

            mfem::DenseMatrix B(6, dof * dim); // In 3D, we have 6 unique strain components.
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

            mfem::DenseMatrix CB(6, dof * dim);
            mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_intpt); // elmat_add is (dof*dim) x (dof*dim)

            elmat.Add(w, elmat_intpt);
        }
    }

    void ThreeDAnisotropicElasticityIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
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

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            mfem::DenseMatrix C;
            stiffness->Eval(C, Trans, ip); // The stiffness matrix is evaluated at each integration point.

            // Here we want to use voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B

            mfem::DenseMatrix B(6, dof * dim); // In 3D, we have 6 unique strain components.
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

            mfem::DenseMatrix CB(6, dof * dim);
            mfem::DenseMatrix elmat_intpt(dof * dim, dof * dim);
            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_intpt); // elmat_add is (dof*dim) x (dof*dim)

            elmat.Add(w, elmat_intpt);
        }
    }

}