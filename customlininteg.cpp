// This file containts custom linear form integrators for the MFEM library..
// List of integrators is given below.

#include "customlininteg.hpp"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void FractureDamageLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;

        shape.SetSize(dof); // vector of size dof
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement
        eldofdamage.SetSize(dof);     // scalar valued damage

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        // Viscosity turned off.
        // for (int i = 0; i < eldofdisp.Size() / dim; i++)
        // {
        //     eldofdamage(i) = (*damage_gf)(eldofs[i]);
        // }
        // Great, now we have damage at each dof. Use that to compute viscosity term at each dof.

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);

        if (ir == NULL)
        {
            ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
        }
        double w, NU, E;

        C.SetSize(str_comp, str_comp);   // Stiffness in Voigt form
        B.SetSize(str_comp, dof * dim);  // Strain displacement matrix
        CB.SetSize(str_comp, dof * dim); // Stiffness times strain displacement
        CBu.SetSize(str_comp);
        Bu.SetSize(str_comp);

        body_pressure.SetSize(str_comp);
        elvect.SetSize(dof);

        C = 0.0;
        B = 0.0;
        body_pressure = 0.0;
        elvect = 0.0;

        double lambda1, lambda2, lambda3;
        double strain_energy(0.0), pressure_energy(0.0), total_energy(0.0);
        double pressure_coeff(0.0);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            if (i == 0)
            {
                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
            } // Constant throughout element

            // Viscosity turned off.
            // eldofdamage *= shape;

            mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            switch (dim)
            {
            case 2:
                // volumetric pressure
                if (volumetric_pressure != nullptr)
                {
                    pressure_coeff = volumetric_pressure->Eval(Tr, ip);
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
                    pressure_coeff = volumetric_pressure->Eval(Tr, ip);
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
            //     total_energy = 0.0;
            // }
            // else
            // {
            //     B.Mult(eldofdisp, Bu); // Bu has dimension strain_comps. This is the strain vector.
            //     strain_energy = mfem::InnerProduct(CBu, Bu);
            // }

            // for now okay, but probably will change it to element average strain energy.
            add(elvect, w * total_energy, shape, elvect);
            // Viscosity turned off.
            // add(elvect, w * viscosity_term, eldofdamage, elvect);
        }
    };
    void FractureHistoryVariableLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {

        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;

        if (elnum == 0)
        {
            shape.SetSize(dof); // vector of size dof
            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            elvect.SetSize(dof);
            elvect = 0.0;

            eldofs.SetSize(dof * dim);    // vector valued for displacement
            eldofdisp.SetSize(dof * dim); // vector valued displacement
            eldofdamage.SetSize(dof);     // scalar valued damage
        }

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*dmg_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
        if (ir == NULL)
        {
            ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
        }
        double w, NU, E;

        if (elnum == 0)
        {
            C.SetSize(str_comp, str_comp);
            B.SetSize(str_comp, dof * dim);
            CB.SetSize(str_comp, dof * dim);
            CBu.SetSize(str_comp);
            Bu.SetSize(str_comp);
        }
        double strain_energy;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            el.CalcPhysShape(Tr, shape);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            // NU = poisson_ratio->Eval(Trans, ip); // need to give coefficients to history variable.
            // E = young_mod->Eval(Trans, ip); // The elastic constants are evaluated at each integration point.

            mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

            // add(elvect, ip.weight * val, shape, elvect);
            if (dim == 2)
            {
                C = 0.0;
                // Plane strain
                // C(0, 0) = C(1, 1) = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU));
                // C(0, 1) = C(1, 0) = E * NU / ((1 + NU) * (1 - 2 * NU));
                // C(2, 2) = E / (2 * (1 + NU));

                // Plane stress
                C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                C(2, 2) = (E * (1 - NU) / (2 * (1 - pow(NU, 2))));

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
            B.Mult(eldofdisp, Bu);   // Bu has dimension strain_comps. This is the strain vector.
            strain_energy = mfem::InnerProduct(CBu, Bu);

            add(elvect, w * strain_energy, shape, elvect); // Hmmm is this all??
        }
    };

    // Eigenstrain body force integrator.
    void EigenstrainBodyForceLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {
        {
            int dof = el.GetDof();
            int dim = el.GetDim();
            int str_comp = (dim == 2) ? 3 : 6;
            int elnum = Tr.ElementNo;

            shape.SetSize(dof); // vector of size dof
            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            elvect.SetSize(dof * dim);
            elvect = 0.0;

            mfem::Vector temp(dof * dim);

            // Great, now we have all components of displacements at each dof.
            // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

            const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
            if (ir == NULL)
            {
                ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
            }

            if (ir == NULL)
            {
                ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
            }
            double w, NU, E;

            C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
            B.SetSize(str_comp, dof * dim); // Strain displacement matrix
            BtC.SetSize(dof * dim, str_comp);
            eps_g.SetSize(str_comp);
            C = 0.0;
            B = 0.0;

            double lambda1, lambda2, lambda3;
            double strain_energy, ip_strain_energy;
            double damage_val, visc_coeff, viscosity_term;

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                el.CalcDShape(ip, dshape);
                Tr.SetIntPoint(&ip);
                el.CalcPhysShape(Tr, shape);
                w = ip.weight * Tr.Weight(); // Quadrature weights

                // damage_val = mfem::InnerProduct(shape, eldofdamage);
                // visc_coeff = viscosity_coeff->Eval(Tr, ip);
                // viscosity_term = damage_val * visc_coeff;

                NU = poisson_ratio->Eval(Tr, ip);
                E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at each integration point.
                growth_strain_coeff->Eval(eps_g, Tr, ip);

                mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use Voigt notation to speed up the assembly process.
                // For this, we need the strain displacement matrix B. The element stiffness can be computed as
                // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
                // The B matrix as 3 rows in 2D and 6 rowd in 3D.

                switch (dim)
                {
                case 2:
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
                        C(2, 2) = (E * (1 - NU) / (2 * (1 - pow(NU, 2))));
                        break;
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
                    break;
                }

                // Now compute the quantity C_{ijkl} u_{k,l} u_{i,j}. Using Voigt notation, of course...
                // This is equivalent to.
                mfem::MultAtB(B, C, BtC); // BtC is (dof * dim) x str_comp
                BtC.Mult(eps_g, temp);

                // for now okay, but probably will change it to element average strain energy.
                add(elvect, w, temp, elvect); // Instead of multiplying the strain_energy, I could add the element average to the vector after.
            }
        }
    }

    void PressureBodyForceLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {
        {
            int dof = el.GetDof();
            int dim = el.GetDim();
            int str_comp = (dim == 2) ? 3 : 6;
            int elnum = Tr.ElementNo;

            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            elvect.SetSize(dof * dim);
            elvect = 0.0;

            mfem::Vector temp(dof * dim);

            // Great, now we have all components of displacements at each dof.
            // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

            const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
            if (ir == NULL)
            {
                ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
            }
            double w, NU, E;

            B.SetSize(str_comp, dof * dim); // Strain displacement matrix
            body_stress.SetSize(str_comp);
            B = 0.0;
            body_stress = 0.0;

            double lambda1, lambda2, lambda3;
            double strain_energy, ip_strain_energy;
            double damage_val, visc_coeff, viscosity_term;
            double pressure_coeff;

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                el.CalcDShape(ip, dshape);
                Tr.SetIntPoint(&ip);
                w = ip.weight * Tr.Weight(); // Quadrature weights

                pressure_coeff = pressure->Eval(Tr, ip);

                mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use Voigt notation to speed up the assembly process.
                // For this, we need the strain displacement matrix B. The element stiffness can be computed as
                // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
                // The B matrix as 3 rows in 2D and 6 rowd in 3D.

                switch (dim)
                {
                case 2:
                    // volumetric pressure
                    body_stress(0) = pressure_coeff;
                    body_stress(1) = pressure_coeff;

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

                    // volumetric pressure
                    body_stress(0) = pressure_coeff;
                    body_stress(1) = pressure_coeff;
                    body_stress(2) = pressure_coeff;

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
                B.MultTranspose(body_stress, temp);
                add(elvect, w, temp, elvect);
            }
        }
    }

    void PressureBodyForceDamageLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {
        {
            int dof = el.GetDof();
            int dim = el.GetDim();
            int str_comp = (dim == 2) ? 3 : 6;
            int elnum = Tr.ElementNo;

            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);
            shape.SetSize(dof);

            eldofs.SetSize(dof); // scalar for damage
            eldofdamage.SetSize(dof);
            eldofs = 0;
            eldofdamage = 0.0;

            elvect.SetSize(dof * dim);
            B.SetSize(str_comp, dof * dim); // Strain displacement matrix
            body_stress.SetSize(str_comp);
            elvec_input.SetSize(dof * dim);
            elvect = 0.0;
            B = 0.0;
            body_stress = 0.0;

            damage_fes->GetElementDofs(elnum, eldofs);

            for (int i = 0; i < eldofdamage.Size(); i++)
            {
                eldofdamage(i) = (*damage_gf)(eldofs[i]);
            }

            // Great, eldofdamage is now set. Now need to use it to construct the interpolated damage value at each quadrature point.
            // i.e., evaluate shape functions at each quadrature point and dot product with eldofdamage to get damage at that point.

            const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
            if (ir == NULL)
            {
                ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
            }
            double w, NU, E;

            double damage_value, degradation_constant;
            double pressure_coeff;

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                el.CalcDShape(ip, dshape);
                Tr.SetIntPoint(&ip);
                el.CalcPhysShape(Tr, shape);

                w = ip.weight * Tr.Weight(); // Quadrature weights

                damage_value = mfem::InnerProduct(shape, eldofdamage);
                degradation_constant = ((1.0 - damage_value) * (1.0 - damage_value)) + k_epsilon; // (1 - d)^{2} + k_{\epsilon}

                pressure_coeff = pressure->Eval(Tr, ip);

                mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use Voigt notation to speed up the assembly process.
                // For this, we need the strain displacement matrix B. The element stiffness can be computed as
                // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
                // The B matrix as 3 rows in 2D and 6 rowd in 3D.

                switch (dim)
                {
                case 2:
                    // volumetric pressure. 2 non-zero components
                    body_stress(0) = pressure_coeff;
                    body_stress(1) = pressure_coeff;

                    // In 2D, we have 3 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        // B(2, spf) = gshape(spf, 1);
                        // B(2, spf + dof) = gshape(spf, 0);
                    }
                    break;
                case 3:

                    // volumetric pressure. 3 non zero components
                    body_stress(0) = pressure_coeff;
                    body_stress(1) = pressure_coeff;
                    body_stress(2) = pressure_coeff;

                    // In 3D, we have 6 unique strain components.
                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf + 2 * dof) = gshape(spf, 2);
                        // B(3, spf + dof) = gshape(spf, 2);
                        // B(3, spf + 2 * dof) = gshape(spf, 1);
                        // B(4, spf) = gshape(spf, 2);
                        // B(4, spf + 2 * dof) = gshape(spf, 0);
                        // B(5, spf) = gshape(spf, 1);
                        // B(5, spf + dof) = gshape(spf, 0);
                    }
                    break;
                }
                B.MultTranspose(body_stress, elvec_input);
                add(elvect, w * degradation_constant, elvec_input, elvect);
            }
        }
    }

    void StVenantKirchoffInternalForceLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
    {
        int dof = el.GetDof();
        int dim = el.GetDim();
        int str_comp = (dim == 2) ? 3 : 6;
        int elnum = Tr.ElementNo;
        int numels = disp_fes->GetNE();

        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        eldofs.SetSize(dof * dim);    // vector valued for displacement
        eldofdisp.SetSize(dof * dim); // vector valued displacement

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);
        if (ir == NULL)
        {
            ir = &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
        }

        double w, NU, E;

        C.SetSize(str_comp, str_comp);           // Stiffness in Voigt form
        BGradDisp.SetSize(dim * dim, dof * dim); // Strain displacement matrix
        Gradu.SetSize(dim * dim);
        Egl.SetSize(str_comp);
        S.SetSize(str_comp);
        F.SetSize(dim * dim);
        BNL.SetSize(str_comp, dof * dim);
        elvec_input.SetSize(dof * dim);
        elvect.SetSize(dof * dim);
        elstrain_ave.SetSize(str_comp);
        elstress_ave.SetSize(str_comp);

        C = 0.0;
        BGradDisp = 0.0;
        BNL = 0.0;
        elvect = 0.0;
        elstrain_ave = 0.0;
        elstress_ave = 0.0;

        int num_int_points = ir->GetNPoints();

        for (int i = 0; i < num_int_points; i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            if (i == 0)
            {
                E = Ey_coeff->Eval(Tr, ip);
                NU = nu_coeff->Eval(Tr, ip); // The elastic constants are evaluated at the first integration point.
            } // Constant throughout element

            mfem::Mult(dshape, Tr.InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use Voigt notation to speed up the assembly process.
            // For this, we need the strain displacement matrix B. The element stiffness can be computed as
            // \int_{\Omega} B^T C B. In Voigt form, the stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            // The B matrix as 3 rows in 2D and 6 rowd in 3D.

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
                if (i == 0)
                {
                    C(0, 0) = C(1, 1) = C(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    C(3, 3) = C(4, 4) = C(5, 5) = E / (2 * (1 + NU));
                }

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
            for (int i = 0; i < dim; i++)
            {
                F(i) += 1.0;
            }

            switch (dim)
            {
            case 3:
                // Fill green lagrange strain vec
                Egl(0) = Gradu(0) + 0.5 * (Gradu(0) * Gradu(0) + Gradu(5) * Gradu(5) + Gradu(7) * Gradu(7));            // E11
                Egl(1) = Gradu(1) + 0.5 * (Gradu(1) * Gradu(1) + Gradu(3) * Gradu(3) + Gradu(8) * Gradu(8));            // E22
                Egl(2) = Gradu(2) + 0.5 * (Gradu(2) * Gradu(2) + Gradu(4) * Gradu(4) + Gradu(6) * Gradu(6));            // E33
                Egl(3) = (Gradu(6) + Gradu(8) + (Gradu(3) * Gradu(4)) + (Gradu(1) * Gradu(6)) + (Gradu(8) * Gradu(2))); // 2 E23
                Egl(4) = (Gradu(4) + Gradu(7) + (Gradu(0) * Gradu(4)) + (Gradu(5) * Gradu(6)) + (Gradu(7) * Gradu(2))); // 2 E13
                Egl(5) = (Gradu(3) + Gradu(5) + (Gradu(0) * Gradu(3)) + (Gradu(5) * Gradu(1)) + (Gradu(7) * Gradu(8))); // 2 E12

                // Strain computed, add to strain_ave.
                elstrain_ave.Add(1.0 / num_int_points, Egl);

                // Compute second piola kirchoff stress
                C.Mult(Egl, S);

                // Stress computed, add to stress_ave.
                elstress_ave.Add(1.0 / num_int_points, S);

                // Fill nonlinear strain displacement matrix B
                // In 3D, we have 6 unique strain components.
                // F: F11, F22, F33, F12, F13, F21, F23, F31, F32
                for (int spf = 0; spf < dof; spf++)
                {
                    // E11
                    BNL(0, spf) = gshape(spf, 0) * F(0);
                    BNL(0, spf + dof) = gshape(spf, 0) * F(5);
                    BNL(0, spf + 2 * dof) = gshape(spf, 0) * F(7);

                    // E22
                    BNL(1, spf) = gshape(spf, 1) * F(3);
                    BNL(1, spf + dof) = gshape(spf, 1) * F(1);
                    BNL(1, spf + 2 * dof) = gshape(spf, 1) * F(8);

                    // E33
                    BNL(2, spf) = gshape(spf, 2) * F(4);
                    BNL(2, spf + dof) = gshape(spf, 2) * F(6);
                    BNL(2, spf + 2 * dof) = gshape(spf, 2) * F(2);

                    // Need to check all of this
                    // E23
                    // BNL(3, spf) = 0.0;
                    BNL(3, spf + dof) = (gshape(spf, 1) * F(6)) + (gshape(spf, 2) * F(1));
                    BNL(3, spf + 2 * dof) = (gshape(spf, 1) * F(2)) + (gshape(spf, 2) * F(8));

                    // E13
                    BNL(4, spf) = (gshape(spf, 0) * F(4)) + (gshape(spf, 2) * F(0));
                    // BNL(4, spf + dof) = 0.0;
                    BNL(4, spf + 2 * dof) = (gshape(spf, 0) * F(2)) + (gshape(spf, 2) * F(7));

                    // E12
                    BNL(5, spf) = (gshape(spf, 0) * F(3)) + (gshape(spf, 1) * F(0));
                    BNL(5, spf + dof) = (gshape(spf, 0) * F(1)) + (gshape(spf, 1) * F(5));
                    // BNL(5, spf + 2 * dof) = 0.0;
                }
                // Assume BNL is assembled correctly. Proceed.
                break;
            }

            // Now compute the quantity S : \bar{E} Using Voigt notation, of course...
            // This is equivalent to.
            BNL.MultTranspose(S, elvec_input);

            add(elvect, -1.0 * w, elvec_input, elvect); // needs to be a negative contribution

            // Add strain_ave and stress_ave to strain_gf and stress_gf in appropriate locations. Element number is known.
            for (int comp = 0; comp < str_comp; comp++)
            {
                (*strain_gf)(elnum + (numels * comp)) = elstrain_ave(comp);
                (*stress_gf)(elnum + (numels * comp)) = elstress_ave(comp);
            }
        }
    }
}
