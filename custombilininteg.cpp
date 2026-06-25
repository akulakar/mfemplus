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

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);

        elmat.SetSize(dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim); // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        C = 0.0;
        B = 0.0;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.
            }

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                if (i == 0)
                {
                    switch (planeApprox)
                    {
                    case 0:
                        // Plane strain
                        C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                        C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;

                    case 1:
                        // Plane stress
                        C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                        C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;
                    }
                }

                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:
                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }

                // In 3D, we have 6 unique strain components.
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
                break;
            }
            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_input); // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w, elmat_input);
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

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);

        elmat.SetSize(dof * dim);
        elmat = 0.0;

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            // int order = Trans.OrderGrad(&el); // correct order?
            // int order = 1;
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        B.SetSize(str_comp, dof * dim); // Strain displacement matrix
        C.SetSize(str_comp, str_comp);  // Elasticity matrix
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        B = 0.0;
        C = 0.0;

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

            stiffness->Eval(C, Trans, ip); // The stiffness matrix is evaluated at each integration point.

            switch (dim)
            {
            case 2:
                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:
                // In 3D, we have 6 unique strain components.
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
                break;
            }

            mfem::Mult(C, B, CB);              // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_input); // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w, elmat_input);
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

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        divshape.SetSize(dim * dof);

        elmat.SetSize(dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        elmat = 0.0;
        B.SetSize(str_comp, dof * dim);     // Strain displacement matrix
        B_vol.SetSize(str_comp, dof * dim); // Strain displacement matrix for volumetric
        B = 0.0;
        B_vol = 0.0;
        elmat_input.SetSize(dof * dim, dof * dim);

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
            B.Add(-1.0 / 3, B_vol);
            mfem::MultAtB(B, B, elmat_input);
            elmat.Add(w * 2.0 * MU, elmat_input);
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
            mfem::QuadratureFunctions1D::GaussLobatto(el.GetOrder() + 2, ir1D);
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

    // Fracture integrators

    void IsotropicElasticityDamageIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        mfem::real_t w, E, NU;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof * dim);
        elmat = 0.0;

        eldofs.SetSize(dof); // scalar for damage
        eldofdamage.SetSize(dof);
        eldofs = 0;
        eldofdamage = 0.0;

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

        C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim); // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim);
        elmat_input.SetSize(dof * dim, dof * dim);
        C = 0.0;
        B = 0.0;
        double damage_value, degradation_constant;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            // k_eps_val not expected to change spatially.
            // multiply degradation constant to the final element matrix.

            el.CalcDShape(ip, dshape);

            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);
            damage_value = mfem::InnerProduct(shape, eldofdamage);                            // Should this be physical or reference shape values? Should be the same.
            degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

            w = ip.weight * Trans.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.
            }
            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                if (i == 0)
                {
                    switch (planeApprox)
                    {
                    case 0:
                        // Plane strain
                        C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                        C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;
                    case 1:
                        // Plane stress
                        C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                        C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;
                    }
                }

                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:

                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }
                // In 3D, we have 6 unique strain components.

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
                break;
            }

            mfem::Mult(C, B, CB);                             // CB is 6 x (dof * dim)
            mfem::MultAtB(B, CB, elmat_input);                // elmat_add is (dof*dim) x (dof*dim)
            elmat.Add(w * degradation_constant, elmat_input); // multiplied by the degradation constant.
        }
    }

    void IsotropicStrainEnergyDamageIntegrator::AssembleElementMatrix(
        const mfem::FiniteElement &el, mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Trans.ElementNo;
        double w, E, NU;

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        shape.SetSize(dof);
        elmat.SetSize(dof);
        elmat = 0.0;

        eldofs.SetSize(dof * dim); // vector valued for displacement
        eldofdisp.SetSize(dof * dim);

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

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);
        body_pressure.SetSize(str_comp);
        elmat_input.SetSize(dof, dof);
        C = 0.0;
        B = 0.0;
        body_pressure = 0.0;

        double lambda1(0.0), lambda2(0.0), lambda3(0.0);
        double strain_energy(0.0), pressure_coeff(0.0), pressure_energy(0.0), total_energy(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Trans.SetIntPoint(&ip);
            el.CalcPhysShape(Trans, shape);

            w = ip.weight * Trans.Weight(); // Quadrature weights

            mfem::Mult(dshape, Trans.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Trans, ip);
                E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at the first integration point.
            }

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                // volumetric pressure. Evaluate at each integration point.
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Trans, ip);
                    body_pressure(0) = pressure_coeff;
                    body_pressure(1) = pressure_coeff;
                }

                if (i == 0)
                {
                    switch (planeApprox)
                    {
                    case 0:
                        // Plane strain
                        C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                        C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;

                    case 1:
                        // Plane stress
                        C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                        C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                        C(2, 2) = E / (2 * (1 + NU));
                        break;
                    }
                }

                // In 2D, we have 3 unique strain components.

                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf) = gshape(spf, 1);
                    B(2, spf + dof) = gshape(spf, 0);
                }
                break;

            case 3:
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Trans, ip);
                    body_pressure(0) = pressure_coeff;
                    body_pressure(1) = pressure_coeff;
                    body_pressure(2) = pressure_coeff;
                }

                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }

                // In 3D, we have 6 unique strain components.
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
                break;
            }

            // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
            // This is equivalent to.
            mfem::Mult(C, B, CB);    // CB is 6 x (dof * dim)
            CB.Mult(eldofdisp, CBu); // CBu has dimension strain_comps. This is the stress vector.
            // For tensile loading, no need for Gershgorin check.

            B.Mult(eldofdisp, Bu);                       // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = mfem::InnerProduct(CBu, Bu); // This is twice the strain energy
            pressure_energy = mfem::InnerProduct(body_pressure, Bu);
            total_energy = strain_energy - 2 * pressure_energy;

            // // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            // if (dim == 2)
            // {
            //     // In 2D lambda min is lambda1.
            //     lambda1 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda2 = (CBu(0) - pressure_coeff + CBu(1) - pressure_coeff) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda3 = lambda1 + 1.0;

            //     // lambda1 = CBu(0) - std::abs(CBu(3));
            //     // lambda2 = CBu(1) - std::abs(CBu(3));
            //     // lambda3 = lambda1 + lambda2; // artificially making it greater than both.
            // }
            // else if (dim == 3)
            // {
            //     lambda1 = CBu(0) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
            //     lambda2 = CBu(1) - pressure_coeff - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
            //     lambda3 = CBu(2) - pressure_coeff - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            // }

            // if (std::min({lambda1, lambda2, lambda3}) < 0)
            // {
            //     total_energy = total_energy * (1e-20);
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = mfem::InnerProduct(CBu, Bu);
            // }

            mfem::AddMult_a_VVt(w * total_energy, shape, elmat); // multiplied by twice the strain energy - twice the pressure energy.
        }
    }

    // ----------------------------------------- Nonlinear elasticity integrators --------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------
    void StVenantKirchoffTangentStiffnessIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;
        mfem::real_t w, E, NU;

        MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }

        elmat.SetSize(dof * dim, dof * dim);

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
        if (ir == NULL)
        {
            int order = 2 * Tr.OrderGrad(&el); // correct order?
            ir = &mfem::IntRules.Get(el.GetGeomType(), order);
        }

        C.SetSize(str_comp, str_comp);                    // Stiffness in Voigt form
        BGradDisp.SetSize(dim * dim, dof * dim);          // Displacement gradient displacement matrix
        BNL.SetSize(str_comp, dof * dim);                 // Nonlinear Displacement Strain matrix
        CB.SetSize(str_comp, dof * dim);                  //
        Sigma.SetSize(dim * dim, dim * dim);              // Matrix form of Second Piola Kirchoff
        elmat_input_1.SetSize(dof * dim, dof * dim);      // element matrix 1
        elmat_input_2.SetSize(dof * dim, dof * dim);      // element matrix 2
        elmat_input_2_temp.SetSize(dof * dim, dof * dim); // element matrix 2 temp
        Gradu.SetSize(dim * dim);                         // Displacement gradient
        F.SetSize(dim * dim);                             // Deformation gradient
        Egl.SetSize(str_comp);                            // Green lagrange strain
        S.SetSize(str_comp);                              // Second Piola Kirchoff

        elmat = 0.0;
        C = 0.0;
        BGradDisp = 0.0;
        BNL = 0.0;
        Sigma = 0.0;
        F = 0.0;
        Gradu = 0.0;
        Egl = 0.0;
        S = 0.0;
        CB = 0.0;
        elmat_input_1 = 0.0;
        elmat_input_2_temp = 0.0;
        elmat_input_2 = 0.0;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);

            Tr.SetIntPoint(&ip);
            w = ip.weight * Tr.Weight();                      // Quadrature weights
            mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at each integration point.
            }

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.
            dim = 3;
            switch (dim)
            {
            case 2:
                if (i == 0)
                {
                    // switch (planeApprox)
                    // {
                    // case 0:
                    // Plane strain
                    C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                    C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                    C(2, 2) = E / (2 * (1 + NU));
                    // break;

                    // case 1:
                    // Plane stress
                    // C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                    // C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                    // C(2, 2) = E / (2 * (1 + NU));
                    // break;
                    // }
                }

                // In 2D, we have 3 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    BGradDisp(0, spf) = gshape(spf, 0);       // u_{1,1}
                    BGradDisp(1, spf + dof) = gshape(spf, 1); // u_{2,2}
                    BGradDisp(2, spf) = gshape(spf, 1);       // u_{1,2}
                    BGradDisp(3, spf + dof) = gshape(spf, 0); // u_{2,1}
                }
                break;

            case 3:
                // if (i == 0)
                // {
                C = 0.0;
                C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                // }

                // In 3D, we have 6 unique strain components.
                for (int spf = 0; spf < dof; spf++)
                {
                    BGradDisp(0, spf) = gshape(spf, 0);           // u_{1,1}
                    BGradDisp(1, spf + dof) = gshape(spf, 1);     // u_{2,2}
                    BGradDisp(2, spf + 2 * dof) = gshape(spf, 2); // u_{3,3}
                    BGradDisp(3, spf) = gshape(spf, 1);           // u_{1,2}
                    BGradDisp(4, spf) = gshape(spf, 2);           // u_{1,3}
                    BGradDisp(5, spf + dof) = gshape(spf, 0);     // u_{2,1}
                    BGradDisp(6, spf + dof) = gshape(spf, 2);     // u_{2,3}
                    BGradDisp(7, spf + 2 * dof) = gshape(spf, 0); // u_{3,1}
                    BGradDisp(8, spf + 2 * dof) = gshape(spf, 1); // u_{3,2}
                }
                break;
            }
            BGradDisp.Mult(eldofdisp, Gradu);
            F = Gradu;
            F = 0.0;
            for (int i = 0; i < dim; i++)
            {
                F(i) += 1.0;
            }

            if (dim == 3)
            {
                // case 3:
                // Fill green lagrange strain vec
                Egl(0) = Gradu(0) + 0.5 * (Gradu(0) * Gradu(0) + Gradu(5) * Gradu(5) + Gradu(7) * Gradu(7));            // E11
                Egl(1) = Gradu(1) + 0.5 * (Gradu(1) * Gradu(1) + Gradu(3) * Gradu(3) + Gradu(8) * Gradu(8));            // E22
                Egl(2) = Gradu(2) + 0.5 * (Gradu(2) * Gradu(2) + Gradu(4) * Gradu(4) + Gradu(6) * Gradu(6));            // E33
                Egl(3) = 0.5 * (Gradu(6) + Gradu(8) + Gradu(3) * Gradu(4) + Gradu(1) * Gradu(6) + Gradu(8) * Gradu(2)); // E23
                Egl(4) = 0.5 * (Gradu(4) + Gradu(7) + Gradu(0) * Gradu(4) + Gradu(5) * Gradu(6) + Gradu(7) * Gradu(2)); // E13
                Egl(5) = 0.5 * (Gradu(3) + Gradu(5) + Gradu(0) * Gradu(3) + Gradu(5) * Gradu(1) + Gradu(7) * Gradu(8)); // E12

                // Compute second piola kirchoff stress
                C.Mult(Egl, S);
                // Matrix representation of second piola for the element matrix.
                Sigma(0, 0) = S(0);
                Sigma(0, 3) = S(5);
                Sigma(0, 4) = S(4);
                Sigma(1, 0) = S(5);
                Sigma(1, 3) = S(1);
                Sigma(1, 4) = S(3);
                Sigma(2, 0) = S(4);
                Sigma(2, 3) = S(3);
                Sigma(2, 4) = S(2);

                Sigma(3, 1) = S(5);
                Sigma(3, 5) = S(0);
                Sigma(3, 6) = S(4);
                Sigma(4, 1) = S(1);
                Sigma(4, 5) = S(5);
                Sigma(4, 6) = S(3);
                Sigma(5, 1) = S(3);
                Sigma(5, 5) = S(4);
                Sigma(5, 6) = S(2);

                Sigma(6, 2) = S(4);
                Sigma(6, 7) = S(0);
                Sigma(6, 8) = S(5);
                Sigma(7, 2) = S(3);
                Sigma(7, 7) = S(5);
                Sigma(7, 8) = S(1);
                Sigma(8, 2) = S(2);
                Sigma(8, 7) = S(4);
                Sigma(8, 8) = S(3);

                // Fill nonlinear strain displacement matrix B
                // In 3D, we have 6 unique strain components.
                // F: F11, F22, F33, F12, F13, F21, F23, F31, F32
                BNL = 0.0;
                // for (int spf = 0; spf < dof; spf++)
                // {
                //     // E11
                //     BNL(0, spf) = gshape(spf, 0) * F(0);
                //     BNL(0, spf + dof) = gshape(spf, 0) * F(5);
                //     BNL(0, spf + 2 * dof) = gshape(spf, 0) * F(7);

                //     // E22
                //     BNL(1, spf) = gshape(spf, 1) * F(3);
                //     BNL(1, spf + dof) = gshape(spf, 1) * F(1);
                //     BNL(1, spf + 2 * dof) = gshape(spf, 1) * F(8);

                //     // E33
                //     BNL(2, spf) = gshape(spf, 2) * F(4);
                //     BNL(2, spf + dof) = gshape(spf, 2) * F(6);
                //     BNL(2, spf + 2 * dof) = gshape(spf, 2) * F(2);

                //     // Need to check all of this
                //     // E23
                //     BNL(3, spf) = 0.0;
                //     BNL(3, spf + dof) = (gshape(spf, 1) * F(6)) + (gshape(spf, 2) * F(1));
                //     BNL(3, spf + 2 * dof) = (gshape(spf, 1) * F(2)) + (gshape(spf, 2) * F(8));

                //     // E13
                //     BNL(4, spf) = (gshape(spf, 0) * F(4)) + (gshape(spf, 2) * F(0));
                //     BNL(4, spf + dof) = 0.0;
                //     BNL(4, spf + 2 * dof) = (gshape(spf, 0) * F(2)) + (gshape(spf, 2) * F(7));

                //     // E12
                //     BNL(5, spf) = (gshape(spf, 0) * F(3)) + (gshape(spf, 1) * F(0));
                //     BNL(5, spf + dof) = (gshape(spf, 0) * F(1)) + (gshape(spf, 1) * F(5));
                //     BNL(5, spf + 2 * dof) = 0.0;
                // }

                BNL = 0.0;
                for (int spf = 0; spf < dof; spf++)
                {
                    BNL(0, spf) = gshape(spf, 0);
                    BNL(1, spf + dof) = gshape(spf, 1);
                    BNL(2, spf + 2 * dof) = gshape(spf, 2);
                    BNL(3, spf + dof) = gshape(spf, 2);
                    BNL(3, spf + 2 * dof) = gshape(spf, 1);
                    BNL(4, spf) = gshape(spf, 2);
                    BNL(4, spf + 2 * dof) = gshape(spf, 0);
                    BNL(5, spf) = gshape(spf, 1);
                    BNL(5, spf + dof) = gshape(spf, 0);
                }
                // Assume BNL is assembled correctly. Proceed.
                // break;
            }

            // elmat input 1
            mfem::MultAtB(BGradDisp, Sigma, elmat_input_1);
            mfem::Mult(elmat_input_1, BGradDisp, elmat_input_1);

            // elmat input 2

            // mfem::MultAtB(BNL, C, elmat_input_2);

            C(0, 0) = C(1, 1) = C(2, 2) = C(3, 3) = C(4, 4) = C(5, 5) = 1.0;
            cout << "BNL norm: " << BNL.FNorm() << endl;
            mfem::Mult(C, BNL, elmat_input_2_temp);
            cout << "elmat 2 temp norm: " << elmat_input_2_temp.FNorm() << endl;
            // mfem::Mult(elmat_input_2, BNL, elmat_input_2);
            // mfem::MultAtB(BNL, elmat_input_2_temp, elmat_input_2);
            mfem::MultAtB(BNL, BNL, elmat_input_2);

            cout << "elmat 2 norm: " << elmat_input_2.FNorm() << endl;

            // elmat = 0.0; // Need to delete.
            // elmat.Add(w, elmat_input_1);
            elmat.Add(w, elmat_input_2);

            // cout << "BNL norm: " << BNL.FNorm() << endl;
            // cout << "Elmat 1 norm:" << elmat_input_1.FNorm() << endl;
            // cout << "Elmat 2 temp norm:" << elmat_input_2_temp.MaxMaxNorm() << endl;
            // cout << "Elmat 2 norm :" << elmat_input_2.MaxMaxNorm() << endl;
            // cout << "Elmat norm:" << elmat.FNorm() << endl;
        }
    }
}
