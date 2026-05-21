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
            if (volumetric_pressure != nullptr)
            {
                pressure_energy = mfem::InnerProduct(body_pressure, Bu);
            }
            total_energy = strain_energy - 2 * pressure_energy;

            // Gershgorin circle theorem for stress. Alternatively, use history variable for strain energy.
            // if (dim == 2)
            // {
            //     // In 2D lambda min is lambda1.
            //     lambda1 = (CBu(0) + CBu(1)) / 2.0 - std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda2 = (CBu(0) + CBu(1)) / 2.0 + std::sqrt(pow((CBu(0) - CBu(1)) / 2.0, 2.0) + pow(CBu(2), 2.0));
            //     lambda3 = lambda1 + 1.0;

            //     // lambda1 = CBu(0) - std::abs(CBu(3));
            //     // lambda2 = CBu(1) - std::abs(CBu(3));
            //     // lambda3 = lambda1 + lambda2; // artificially making it greater than both.
            // }
            // else if (dim == 3)
            // {
            //     lambda1 = CBu(0) - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
            //     lambda2 = CBu(1) - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
            //     lambda3 = CBu(2) - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            // }

            // if (std::min({lambda1, lambda2, lambda3}) > 0)
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = mfem::InnerProduct(CBu, Bu);
            // }
            // else
            //     strain_energy = 1.0e-20;

            mfem::AddMult_a_VVt(w * total_energy, shape, elmat); // multiplied by twice the strain energy - twice the pressure energy.
        }
    }
}
