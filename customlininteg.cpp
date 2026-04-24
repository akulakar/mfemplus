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

        shape.SetSize(dof); // vector of size dof
        dshape.SetSize(dof, dim);
        gshape.SetSize(dof, dim);
        elvect.SetSize(dof);
        elvect = 0.0;

        mfem::Array<int> eldofs(dof * dim); // vector valued for displacement
        mfem::Vector eldofdisp(dof * dim);

        int elnum = Tr.ElementNo;

        disp_fes->GetElementVDofs(elnum, eldofs);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            eldofdisp(i) = (*disp_gf)(eldofs[i]);
        }
        // Great, now we have all components of displacements at each dof.
        // Next, construct the stiffness matrix C, compute displacement gradients, and take inner product.

        const mfem::IntegrationRule *ir = GetIntegrationRule(el, Tr);

        if (ir == NULL)
        {
            // ir = &IntRules.Get(el.GetGeomType(),
            //                    oa * el.GetOrder() + ob + Tr.OrderW());
            ir = &mfem::IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);
        }
        double w, NU, E;

        mfem::DenseMatrix C(str_comp, str_comp);   // Stiffness in Voigt form
        mfem::DenseMatrix B(str_comp, dof * dim);  // Strain displacement matrix
        mfem::DenseMatrix CB(str_comp, dof * dim); // Stiffness times strain displacement
        mfem::Vector CBu(str_comp);
        mfem::Vector Bu(str_comp);
        double lambda1, lambda2, lambda3;
        double strain_energy;

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el.CalcDShape(ip, dshape);
            Tr.SetIntPoint(&ip);
            w = ip.weight * Tr.Weight(); // Quadrature weights

            NU = poisson_ratio->Eval(Tr, ip);
            E = young_mod->Eval(Tr, ip); // The elastic constants are evaluated at each integration point.

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
            B.Mult(eldofdisp, Bu);   // Bu has dimension strain_comps. This is the strain vector.

            if (dim == 3)
            {
                lambda1 = CBu(0) - std::abs(CBu(5)) - std::abs(CBu(4)); // \lambda_{1} = \sigma_{11} - |\sigma_{12}| - |\sigma_{13}|
                lambda2 = CBu(1) - std::abs(CBu(5)) - std::abs(CBu(3)); // \lambda_{2} = \sigma_{22} - |\sigma_{12}| - |\sigma_{23}|
                lambda3 = CBu(2) - std::abs(CBu(4)) - std::abs(CBu(3)); // \lambda_{3} = \sigma_{33} - |\sigma_{13}| - |\sigma_{23}|
            }

            if (std::min({lambda1, lambda2, lambda3}) > 0)
            {
                strain_energy = mfem::InnerProduct(CBu, Bu);
            }
            else
                strain_energy = 0.0;

            el.CalcPhysShape(Tr, shape);
            add(elvect, w * strain_energy, shape, elvect);
        }
    };
}
