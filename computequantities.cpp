// This constains the implementation of custom functions written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
// List of functions are given below.
// 1. Elasticity
//      (a) Recover strain tensor for each element using B matrix given a displacement grid function.
//      (b) Recover stress tensor for each element given strain and constitutive law.
//      (c) Recover divergence of displacement field using B matrix.
//      (d) Recover curl of displacement field using B matrix.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#include "computequantities.hpp"
// #include "omp.h"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void ElementStressStrain::ComputeElementStrain(mfem::GridFunction &disp,
                                                   int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstrain)
    {
        // To recover element strain, we need the B matrix at each integration point multiplied by
        // the displacement vector at each node in the element. Therefore, the strain varies
        // within the element if the B matrix is not constant. However, the strain is averaged
        // within an element.

        AccessMFEMFunctions accessfunc;

        const mfem::FiniteElement *element = fes->GetFE(elnum);
        mfem::ElementTransformation *Trans = fes->GetElementTransformation(elnum);

        int dof = element->GetDof();
        int dim = element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Array<double> eldofdisp(dof * dim);

        fes->GetElementVDofs(elnum, eldofs);
        int eltype = fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp[i] = disp(dof);
        }

        mfem::real_t w, E, NU;
        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);

        elstrain.SetSize(dim == 2 ? 3 : 6);

        elstrain = 0.0;

        const mfem::IntegrationRule *ir = accessfunc.GetIntegrationRule(*element, *Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans->OrderGrad(element); // correct order?
            ir = &mfem::IntRules.Get(element->GetGeomType(), order);
        }

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.

            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Here we want to use voigt notation to compute the strain vector.
            // The strain vector has 3 components in 2D and 6 components in 3D.
            mfem::DenseMatrix B;

            if (dim == 2)
            {
                B.SetSize(3, dof * dim); // In 2D, we have 3 unique strain components.
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
                B.SetSize(6, dof * dim); // In 3D, we have 6 unique strain components.
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

            mfem::Vector temp(elstrain.Size());
            temp = 0.0;
            B.Mult(eldofdisp, temp);
            double w = 1.0 / (ir->GetNPoints());
            elstrain.Add(w, temp);
        }
    };

    void ElementStressStrain::ComputeElementStress(mfem::Vector &elstrain, mfem::Coefficient &e, mfem::Coefficient &nu,
                                                   int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress)
    {
        AccessMFEMFunctions accessfunc;

        const mfem::FiniteElement *element = fes->GetFE(elnum);
        mfem::ElementTransformation *Trans = fes->GetElementTransformation(elnum);
        int dim = element->GetDim();
        // Stiffness in Voigt notation. The stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
        mfem::DenseMatrix C;

        if (dim == 2)
        {
            C.SetSize(3, 3);
            elstress.SetSize(3);
        }
        else if (dim == 3)
        {
            C.SetSize(6, 6);
            elstress.SetSize(6);
        }

        C = 0.0;
        elstress = 0.0;

        const mfem::IntegrationRule *ir = accessfunc.GetIntegrationRule(*element, *Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans->OrderGrad(element); // correct order?
            ir = &mfem::IntRules.Get(element->GetGeomType(), order);
        }

        mfem::real_t NU, E;
        mfem::DenseMatrix Ctemp;
        Ctemp.SetSize(C.NumRows(), C.NumCols());
        double w = 1.0 / ir->GetNPoints();

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            NU = E = 0;
            Ctemp = 0.0;
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            Trans->SetIntPoint(&ip);

            // The elastic constants are evaluated at each integration point.
            mfem::real_t NU_temp, E_temp;
            NU = nu.Eval(*Trans, ip);
            E = e.Eval(*Trans, ip);

            if (dim == 2)
            {
                // Plane strain

                C(0, 0) = C(1, 1) = E / (1 - pow(NU, 2));
                C(0, 1) = C(1, 0) = (E * NU) / (1 - pow(NU, 2));
                C(2, 2) = (E * (1 - NU)) / (2 * (1 - pow(NU, 2)));

                // Plane stress

                // C(0, 0) = C(1, 1) = (E / (1 - pow(NU, 2)));
                // C(0, 1) = C(1, 0) = (E * NU / (1 - pow(NU, 2)));
                // C(2, 2) = (E * (1 - NU) / (2 * (1 - pow(NU, 2))));

                // 3D to 2D, but this is improper

                // C(0, 0) = C(1, 1) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                // C(0, 1) = C(1, 0) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                // C(2, 2) =  E / (2 * (1 + NU));
            }
            else if (dim == 3)
            {
                Ctemp(0, 0) = Ctemp(1, 1) = Ctemp(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                Ctemp(0, 1) = Ctemp(0, 2) = Ctemp(1, 0) = Ctemp(1, 2) = Ctemp(2, 0) = Ctemp(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                Ctemp(3, 3) = Ctemp(4, 4) = Ctemp(5, 5) = E / (2 * (1 + NU));
            }

            C.Add(w, Ctemp);
        };

        C.Mult(elstrain, elstress);
    };

    void ElementStressStrain::ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat,
                                                   int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress)
    {
        AccessMFEMFunctions accessfunc;

        const mfem::FiniteElement *element = fes->GetFE(elnum);
        mfem::ElementTransformation *Trans = fes->GetElementTransformation(elnum);
        int dim = element->GetDim();
        // Stiffness in Voigt notation. The stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
        mfem::DenseMatrix C;

        if (dim == 2)
            C.SetSize(3, 3);
        else if (dim == 3)
            C.SetSize(6, 6);

        const mfem::IntegrationRule *ir = accessfunc.GetIntegrationRule(*element, *Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans->OrderGrad(element); // correct order?
            ir = &mfem::IntRules.Get(element->GetGeomType(), order);
        }

        double w = 1.0 / ir->GetNPoints();
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            Trans->SetIntPoint(&ip);

            mfem::DenseMatrix C_temp;
            Cmat.Eval(C_temp, *Trans, ip);

            C.Add(w, C_temp);
        }
        C.Mult(elstrain, elstress);
    };

    void ElementStressStrain::ComputeBoundaryElementArea(mfem::real_t &area, const mfem::FiniteElement *element,
                                                         mfem::ElementTransformation *Trans)
    {
        area = 0.0;
        mfem::real_t w = 0.0;
        int order = element->GetOrder();
        int geometry = element->GetGeomType();

        const mfem::IntegrationRule *ir = &mfem::IntRules.Get(geometry, order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            Trans->SetIntPoint(&ip);
            w = ip.weight * Trans->Weight(); // Quadrature weights
            area += w;
        }
    };

    void ElementStressStrain::ComputeElementDilatation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::real_t &dilatation)
    {
        AccessMFEMFunctions accessfunc;
        // The dilation is div(u). We need the B matrix at each integration point multiplied by
        // the displacement vector at each node in the element. Therefore, the dilatation varies
        // within the element if the B matrix is not constant. However, the dilatation is averaged
        // within an element.

        const mfem::FiniteElement *element = fes->GetFE(elnum);
        mfem::ElementTransformation *Trans = fes->GetElementTransformation(elnum);

        int dof = element->GetDof();
        int dim = element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Array<double> eldofdisp(dof * dim);

        fes->GetElementVDofs(elnum, eldofs);
        int eltype = fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp[i] = (*disp)(dof);
        }

        mfem::real_t w;
        mfem::real_t w_sum = 0.0;
        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);

        dilatation = 0.0;

        const mfem::IntegrationRule *ir = accessfunc.GetIntegrationRule(*element, *Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans->OrderGrad(element); // correct order?
            ir = &mfem::IntRules.Get(element->GetGeomType(), order);
        }

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            w = ip.weight * Trans->Weight();
            w_sum += w;
            element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.

            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Arrange the gradients of the shape functions in the B matrix to recover the scalar div(u).
            mfem::DenseMatrix B;

            if (dim == 2)
            {
                B.SetSize(2, dof * dim); // In 2D, we have 2 dilatational components.
                B = 0.0;
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                }
            }

            else if (dim == 3)
            {
                B.SetSize(3, dof * dim); // In 3D, we have 3 dilatational components.
                B = 0.0;
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf) = gshape(spf, 0);
                    B(1, spf + dof) = gshape(spf, 1);
                    B(2, spf + 2 * dof) = gshape(spf, 2);
                }
            }
            mfem::Vector temp(B.NumRows());
            mfem::real_t dilatation_temp = 0.0;
            temp = 0.0;
            B.Mult(eldofdisp, temp);
            dilatation_temp = temp.Sum();
            dilatation += w * dilatation_temp;
        }
        dilatation *= 1.0 / w_sum;
    };

    void ElementStressStrain::ComputeElementRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &rotation)
    {
        AccessMFEMFunctions accessfunc;
        // The  is div(u). We need the B matrix at each integration point multiplied by
        // the displacement vector at each node in the element. Therefore, the dilatation varies
        // within the element if the B matrix is not constant. However, the dilatation is averaged
        // within an element.

        const mfem::FiniteElement *element = fes->GetFE(elnum);
        mfem::ElementTransformation *Trans = fes->GetElementTransformation(elnum);

        int dof = element->GetDof();
        int dim = element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Array<double> eldofdisp(dof * dim);

        fes->GetElementVDofs(elnum, eldofs);
        int eltype = fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp[i] = (*disp)(dof);
        }

        mfem::real_t w;
        mfem::real_t w_sum = 0.0;
        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);

        rotation.SetSize(dim == 2 ? 1 : 3);
        rotation = 0.0;

        const mfem::IntegrationRule *ir = accessfunc.GetIntegrationRule(*element, *Trans);
        if (ir == NULL)
        {
            int order = 2 * Trans->OrderGrad(element); // correct order?
            ir = &mfem::IntRules.Get(element->GetGeomType(), order);
        }

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            w = ip.weight * Trans->Weight();
            w_sum += w;
            element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.

            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            // Arrange the gradients of the shape functions in the B matrix to recover the scalar rot(u) in 2D and the vector curl(u) in 3D.
            mfem::DenseMatrix B;

            if (dim == 2)
            {
                B.SetSize(1, dof * dim); // In 2D, we have a rot scalar.
                B = 0.0;
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf + dof) = gshape(spf, 0);
                    B(0, spf) = -gshape(spf, 1);
                }
            }

            else if (dim == 3)
            {
                B.SetSize(3, dof * dim); // In 3D, we have a curl vector with 3 components.
                B = 0.0;
                for (int spf = 0; spf < dof; spf++)
                {
                    B(0, spf + 2 * dof) = gshape(spf, 1);
                    B(0, spf + dof) = -gshape(spf, 2);
                    B(1, spf + 2 * dof) = -gshape(spf, 0);
                    B(1, spf) = gshape(spf, 2);
                    B(2, spf + dof) = gshape(spf, 0);
                    B(2, spf) = -gshape(spf, 1);
                }
            }
            mfem::Vector temp(B.NumRows());
            B.Mult(eldofdisp, temp);
            add(rotation, w, temp, rotation);
        }
        rotation *= (1.0 / 2.0) * (1.0 / w_sum);
    }

    //------------------------------------------------------------------------------------------------------------------------------------
    // General function to compute displacement gradients at each integration point.
    mfem::Vector ElementStressStrain::ComputeElementDisplacementGradients(mfem::DenseMatrix &gshape, mfem::Vector &eldofdisp)
    {
        // Arrange the gradients of the shape functions in the B matrix to recover the displacement gradients.
        // There are 9 displacement gradients in 3D and 4 in 2D.

        int dim = gshape.NumCols();
        int dof = gshape.NumRows();
        mfem::DenseMatrix B;
        mfem::Vector temp;

        if (dim == 2)
        {
            B.SetSize(4, dof * dim);
            temp.SetSize(4);
            B = 0.0;
            for (int spf = 0; spf < dof; spf++)
            {
                B(0, spf) = gshape(spf, 0);       // u_{1,1}
                B(1, spf + dof) = gshape(spf, 1); // u_{2,2}
                B(2, spf) = gshape(spf, 1);       // u_{1,2}
                B(3, spf + dof) = gshape(spf, 0); // u_{2,1}
            }
        }

        else if (dim == 3)
        {
            B.SetSize(9, dof * dim);
            temp.SetSize(9);
            B = 0.0;
            for (int spf = 0; spf < dof; spf++)
            {
                B(0, spf) = gshape(spf, 0);           // u_{1,1}
                B(1, spf + dof) = gshape(spf, 1);     // u_{2,2}
                B(2, spf + 2 * dof) = gshape(spf, 2); // u_{3,3}
                B(3, spf) = gshape(spf, 1);           // u_{1,2}
                B(4, spf) = gshape(spf, 2);           // u_{1,3}
                B(5, spf + dof) = gshape(spf, 0);     // u_{2,1}
                B(6, spf + dof) = gshape(spf, 2);     // u_{2,3}
                B(7, spf + 2 * dof) = gshape(spf, 0); // u_{3,1}
                B(8, spf + 2 * dof) = gshape(spf, 1); // u_{3,2}
            }
        }
        B.Mult(eldofdisp, temp);
        return temp;
    };

    void ElementStressStrain::ComputeElementStrainip(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &elstrain)
    {
        AccessMFEMFunctions accessfunc;
        const mfem::FiniteElement *disp_element = disp_fes->GetFE(elnum);
        const mfem::FiniteElement *L2_element = L2_fes->GetFE(elnum);

        mfem::ElementTransformation *Trans = disp_fes->GetElementTransformation(elnum);

        int dof = disp_element->GetDof();
        int dim = disp_element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Vector eldofdisp(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);
        int eltype = disp_fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int loc_dof = eldofs[i];
            eldofdisp(i) = disp(loc_dof);
        }

        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);

        mfem::Vector ip_strain, ip_rotation;
        ip_strain.SetSize(dim == 2 ? 3 : 6);
        ip_rotation.SetSize(dim == 2 ? 1 : 3);
        ip_strain = 0.0;
        ip_rotation = 0.0;

        const mfem::IntegrationRule *ir(&(L2_element->GetNodes()));
        int num_int_points = ir->GetNPoints();

        elstrain.SetSize(dim == 2 ? (3 * num_int_points) : (6 * num_int_points));
        elstrain = 0.0;

        for (int i = 0; i < num_int_points; i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            disp_element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.
            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            mfem::Vector disp_gradients(ComputeElementDisplacementGradients(gshape, eldofdisp)); // Compute displacement gradients vector.
            // Now, sum up the displacement gradients properly to get strain and curl(u) (or rot(u) in 2D).
            if (dim == 2)
            {
                elstrain(i) = disp_gradients(0);                                                  // u_{1,1}
                elstrain(i + num_int_points) = disp_gradients(1);                                 // u_{2,2}
                elstrain(i + 2 * num_int_points) = 0.5 * (disp_gradients(2) + disp_gradients(3)); // \frac{1}{2} (u_{1,2} + u_{2,1})
            }

            if (dim == 3)
            {
                elstrain(i) = disp_gradients(0);                                                  // u_{1,1}
                elstrain(i + num_int_points) = disp_gradients(1);                                 // u_{2,2}
                elstrain(i + 2 * num_int_points) = disp_gradients(2);                             // u_{3,3}
                elstrain(i + 3 * num_int_points) = 0.5 * (disp_gradients(6) + disp_gradients(8)); // \frac{1}{2} (u_{2,3} + u_{3,2})
                elstrain(i + 4 * num_int_points) = 0.5 * (disp_gradients(4) + disp_gradients(7)); // \frac{1}{2} (u_{1,3} + u_{3,1})
                elstrain(i + 5 * num_int_points) = 0.5 * (disp_gradients(3) + disp_gradients(5)); // \frac{1}{2} (u_{1,2} + u_{2,1})
            }
        }
    };

    // Function to compute rotation and strain. However, not averaged, but each integration point is a dof in an L2GridFunction.
    void ElementStressStrain::ComputeElementStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &strain, mfem::Vector &rotation)
    {
        AccessMFEMFunctions accessfunc;
        const mfem::FiniteElement *disp_element = disp_fes->GetFE(elnum);
        const mfem::FiniteElement *L2_element = L2_fes->GetFE(elnum);

        mfem::ElementTransformation *Trans = disp_fes->GetElementTransformation(elnum);

        int dof = disp_element->GetDof();
        int dim = disp_element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Vector eldofdisp(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);
        int eltype = disp_fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp(i) = (*disp)(dof);
        }

        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);

        mfem::Vector ip_strain, ip_rotation;
        ip_strain.SetSize(dim == 2 ? 3 : 6);
        ip_rotation.SetSize(dim == 2 ? 1 : 3);
        ip_strain = 0.0;
        ip_rotation = 0.0;

        const mfem::IntegrationRule *ir(&(L2_element->GetNodes()));
        int num_int_points = ir->GetNPoints();

        strain.SetSize(dim == 2 ? (3 * num_int_points) : (6 * num_int_points));
        rotation.SetSize(dim == 2 ? (1 * num_int_points) : (3 * num_int_points));
        strain = 0.0;
        rotation = 0.0;

        for (int i = 0; i < num_int_points; i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            disp_element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.
            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            mfem::Vector disp_gradients(ComputeElementDisplacementGradients(gshape, eldofdisp)); // Compute displacement gradients vector.
            // Now, sum up the displacement gradients properly to get strain and curl(u) (or rot(u) in 2D).
            if (dim == 2)
            {
                strain(i) = disp_gradients(0);                                                  // u_{1,1}
                strain(i + num_int_points) = disp_gradients(1);                                 // u_{2,2}
                strain(i + 2 * num_int_points) = 0.5 * (disp_gradients(2) + disp_gradients(3)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                rotation(i) = disp_gradients(3) - disp_gradients(2); // u_{2,1} - u_{1,2}
            }

            if (dim == 3)
            {
                strain(i) = disp_gradients(0);                                                  // u_{1,1}
                strain(i + num_int_points) = disp_gradients(1);                                 // u_{2,2}
                strain(i + 2 * num_int_points) = disp_gradients(2);                             // u_{3,3}
                strain(i + 3 * num_int_points) = 0.5 * (disp_gradients(6) + disp_gradients(8)); // \frac{1}{2} (u_{2,3} + u_{3,2})
                strain(i + 4 * num_int_points) = 0.5 * (disp_gradients(4) + disp_gradients(7)); // \frac{1}{2} (u_{1,3} + u_{3,1})
                strain(i + 5 * num_int_points) = 0.5 * (disp_gradients(3) + disp_gradients(5)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                rotation(i) = 0.5 * (disp_gradients(8) - disp_gradients(6));                      // \frac{1}{2} u_{3,2} - u_{2,3}
                rotation(i + num_int_points) = 0.5 * (disp_gradients(4) - disp_gradients(7));     // \frac{1}{2} - u_{3,1} + u_{1,3}
                rotation(i + 2 * num_int_points) = 0.5 * (disp_gradients(5) - disp_gradients(3)); // \frac{1}{2} u_{2,1} - u_{1,2}
            }
        }
    }

    void ElementStressStrain::ComputeElementMaxShearStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &max_shear_strain, mfem::Vector &rotation)
    {
        AccessMFEMFunctions accessfunc;
        const mfem::FiniteElement *disp_element = disp_fes->GetFE(elnum);
        const mfem::FiniteElement *L2_element = L2_fes->GetFE(elnum);

        mfem::ElementTransformation *Trans = disp_fes->GetElementTransformation(elnum);

        int dof = disp_element->GetDof();
        int dim = disp_element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Vector eldofdisp(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);
        int eltype = disp_fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp(i) = (*disp)(dof);
        }

        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);
        const mfem::IntegrationRule *ir(&(L2_element->GetNodes()));
        int num_int_points = ir->GetNPoints();

        max_shear_strain.SetSize(num_int_points);
        rotation.SetSize(dim == 2 ? (num_int_points) : (3 * num_int_points));
        max_shear_strain = 0.0;
        rotation = 0.0;

        for (int i = 0; i < num_int_points; i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            disp_element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.
            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            mfem::Vector disp_gradients(ComputeElementDisplacementGradients(gshape, eldofdisp)); // Compute displacement gradients vector.
            // Now, sum up the displacement gradients properly to get strain and curl(u) (or rot(u) in 2D).
            if (dim == 2)
            {
                double eps11 = disp_gradients(0);                             // u_{1,1}
                double eps22 = disp_gradients(1);                             // u_{2,2}
                double eps12 = 0.5 * (disp_gradients(2) + disp_gradients(3)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                double max_strain_val = pow(pow((eps11 - eps22) / 2.0, 2) + pow(eps12, 2), 0.5);

                max_shear_strain(i) = abs(max_strain_val);

                rotation(i) = disp_gradients(3) - disp_gradients(2); // u_{2,1} - u_{1,2}
            }

            if (dim == 3)
            {
                // This is certainly wrong. Needs to be changed.
                double eps11 = disp_gradients(0);                             // u_{1,1}
                double eps22 = disp_gradients(1);                             // u_{2,2}
                double eps12 = 0.5 * (disp_gradients(2) + disp_gradients(3)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                double max_strain_val = +pow(pow((eps11 - eps22) / 2.0, 2) + pow(eps12, 2), 0.5);
                double min_strain_val = -pow(pow((eps11 - eps22) / 2.0, 2) + pow(eps12, 2), 0.5);

                max_shear_strain(i) = (abs(max_strain_val) > abs(min_strain_val)) ? max_strain_val : min_strain_val;
                rotation(i) = 0.5 * (disp_gradients(8) - disp_gradients(6));                      // \frac{1}{2} u_{3,2} - u_{2,3}
                rotation(i + num_int_points) = 0.5 * (disp_gradients(4) - disp_gradients(7));     // \frac{1}{2} - u_{3,1} + u_{1,3}
                rotation(i + 2 * num_int_points) = 0.5 * (disp_gradients(5) - disp_gradients(3)); // \frac{1}{2} u_{2,1} - u_{1,2}
            }
        }
    }

    void ElementStressStrain::ComputeElementStrainMaxShearStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &strain, mfem::Vector &max_shear_strain, mfem::Vector &rotation)
    {
        AccessMFEMFunctions accessfunc;
        const mfem::FiniteElement *disp_element = disp_fes->GetFE(elnum);
        const mfem::FiniteElement *L2_element = L2_fes->GetFE(elnum);

        mfem::ElementTransformation *Trans = disp_fes->GetElementTransformation(elnum);

        int dof = disp_element->GetDof();
        int dim = disp_element->GetDim();

        mfem::Array<int> eldofs(dof * dim);
        mfem::Vector eldofdisp(dof * dim);

        disp_fes->GetElementVDofs(elnum, eldofs);
        int eltype = disp_fes->GetElementType(elnum);

        for (int i = 0; i < eldofdisp.Size(); i++)
        {
            int dof = eldofs[i];
            eldofdisp(i) = (*disp)(dof);
        }

        mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);
        const mfem::IntegrationRule *ir(&(L2_element->GetNodes()));
        int num_int_points = ir->GetNPoints();

        strain.SetSize(dim == 2 ? (num_int_points * 3) : (num_int_points * 6));
        max_shear_strain.SetSize(num_int_points);
        rotation.SetSize(dim == 2 ? (num_int_points) : (3 * num_int_points));
        strain = 0.0;
        max_shear_strain = 0.0;
        rotation = 0.0;

        for (int i = 0; i < num_int_points; i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);
            disp_element->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.
            Trans->SetIntPoint(&ip);
            mfem::Mult(dshape, Trans->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

            mfem::Vector disp_gradients(ComputeElementDisplacementGradients(gshape, eldofdisp)); // Compute displacement gradients vector.
            // Now, sum up the displacement gradients properly to get strain and curl(u) (or rot(u) in 2D).
            if (dim == 2)
            {
                double eps11 = disp_gradients(0);                             // u_{1,1}
                double eps22 = disp_gradients(1);                             // u_{2,2}
                double eps12 = 0.5 * (disp_gradients(2) + disp_gradients(3)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                strain(3 * i) = eps11;
                strain(3 * i + 1) = eps22;
                strain(3 * i + 2) = eps12;

                double max_strain_val = pow(pow((eps11 - eps22) / 2.0, 2) + pow(eps12, 2), 0.5);

                max_shear_strain(i) = abs(max_strain_val);

                rotation(i) = disp_gradients(3) - disp_gradients(2); // u_{2,1} - u_{1,2}
            }

            if (dim == 3)
            {
                double eps11 = disp_gradients(0);                             // u_{1,1}
                double eps22 = disp_gradients(1);                             // u_{2,2}
                double eps33 = disp_gradients(2);                             // u_{3,3}
                double eps23 = 0.5 * (disp_gradients(6) + disp_gradients(8)); // \frac{1}{2} (u_{2,3} + u_{3,2}
                double eps13 = 0.5 * (disp_gradients(5) + disp_gradients(7)); // \frac{1}{2} (u_{1,3} + u_{3,1}
                double eps12 = 0.5 * (disp_gradients(3) + disp_gradients(5)); // \frac{1}{2} (u_{1,2} + u_{2,1})

                strain(3 * i) = eps11;
                strain(3 * i + 1) = eps22;
                strain(3 * i + 2) = eps33;
                strain(3 * i + 3) = eps23;
                strain(3 * i + 4) = eps13;
                strain(3 * i + 5) = eps12;

                // This is certainly wrong. Needs to be changed.
                double max_shear_strain_val = pow(pow((eps11 - eps22) / 2.0, 2) + pow(eps12, 2), 0.5);
                max_shear_strain(i) = max_shear_strain_val;

                rotation(i) = 0.5 * (disp_gradients(8) - disp_gradients(6));                      // \frac{1}{2} u_{3,2} - u_{2,3}
                rotation(i + num_int_points) = 0.5 * (disp_gradients(4) - disp_gradients(7));     // \frac{1}{2} - u_{3,1} + u_{1,3}
                rotation(i + 2 * num_int_points) = 0.5 * (disp_gradients(5) - disp_gradients(3)); // \frac{1}{2} u_{2,1} - u_{1,2}
            }
        }
    };
    //------------------------------------------------------------------------------------------------------------------------------------

    void GlobalStressStrain::GlobalStrain(mfem::GridFunction &disp, mfem::GridFunction &strain)
    {
        int numels = disp_fespace->GetNE();
        int dim = mesh->Dimension();
        // There are 3 strain components in 2D and 6 in 3D.
        int str_comp = (dim == 2) ? 3 : 6;

        // Now start a loop that assembles strain vector element by element, and assembles the grid function.
        // The zeroth to numels index is \epsilon_{11}, then \epsilon_{22}, \epsilon_{33},
        // \epsilon_{23}, \epsilon_{13}, \epsilon_{12}.

        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector elstrain;
            ElementStressStrain Element;
            Element.ComputeElementStrain(disp, elnum, disp_fespace, elstrain);
            for (int comp = 0; comp < str_comp; comp++)
                strain(elnum + (numels * comp)) = elstrain(comp);
        }
    };

    void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e,
                                          mfem::Coefficient &nu, mfem::GridFunction &stress)
    {

        int numels = disp_fespace->GetNE();
        int dim = mesh->Dimension();
        // There are 3 strain components in 2D and 6 in 3D.
        int str_comp = (dim == 2) ? 3 : (dim == 3) ? 6
                                                   : 3;

        stress.SetSize(strain.Size());
        stress = 0.0;

        // Now start a loop that assembles stress vector element by element, and assembles the grid function.
        // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33},
        // \sigma_{23}, \sigma_{13}, \sigma_{12}.

        for (int elnum = 0; elnum < numels; elnum++)
        {

            mfem::Vector elstress(str_comp);
            mfem::Vector elstrain(str_comp);

            for (int comp = 0; comp < str_comp; comp++)
                elstrain(comp) = strain(elnum + (numels * comp));

            ElementStressStrain Element;
            Element.ComputeElementStress(elstrain, e, nu, elnum, disp_fespace, elstress);

            for (int comp = 0; comp < str_comp; comp++)
                stress(elnum + (numels * comp)) = elstress(comp);
        }
    };

    void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain,
                                          mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress)
    {

        int numels = disp_fespace->GetNE();
        int dim = mesh->Dimension();
        // There are 3 strain components in 2D and 6 in 3D.
        int str_comp = (dim == 2) ? 3 : (dim == 3) ? 6
                                                   : 3;

        stress.SetSize(strain.Size());
        stress = 0.0;

        // Now start a loop that assembles stress vector element by element, and assembles the grid function.
        // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33},
        // \sigma_{23}, \sigma_{13}, \sigma_{12}.

        for (int elnum = 0; elnum < numels; elnum++)
        {

            mfem::Vector elstress(str_comp);
            mfem::Vector elstrain(str_comp);

            for (int comp = 0; comp < str_comp; comp++)
                elstrain(comp) = strain(elnum + (numels * comp));

            ElementStressStrain Element;
            Element.ComputeElementStress(elstrain, Cmat, elnum, disp_fespace, elstress);

            for (int comp = 0; comp < str_comp; comp++)
                stress(elnum + (numels * comp)) = elstress(comp);
        }
    };

    double GlobalStressStrain::ComputeBoundaryForce(mfem::GridFunction &stress, int &bdr_attribute, int &component)
    {
        double top_force = 0.0;

        int num_bdr_els = disp_fespace->GetNBE();
        int num_els = disp_fespace->GetNE();
        for (int i = 0; i < num_bdr_els; i++)
        {
            int bdr_attr = disp_fespace->GetBdrAttribute(i);
            if (bdr_attr == bdr_attribute)
            {
                mfem::FaceElementTransformations *ftr = disp_fespace->GetMesh()->GetBdrFaceTransformations(i);
                int volume_elem_index = ftr->Elem1No;
                double element_sig = stress(volume_elem_index + (component - 1) * num_els);
                double element_force = 0.0;
                mfem::real_t element_area;
                const mfem::FiniteElement *surf_el = disp_fespace->GetBE(i);
                mfem::ElementTransformation *el_trans = disp_fespace->GetBdrElementTransformation(i);
                mfemplus::ElementStressStrain ElementArea;
                ElementArea.ComputeBoundaryElementArea(element_area, surf_el, el_trans);
                element_force = element_area * element_sig;
                top_force += element_force;
            }
        }
        return top_force;
    };

    void GlobalStressStrain::GlobalDilatation(mfem::GridFunction *disp, mfem::GridFunction *dilatation)
    {
        int numels = disp_fespace->GetNE();
        int dim = mesh->Dimension();
        // There is one measure of dilatation per element.
        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::real_t element_dilatation;
            ElementStressStrain Element;
            Element.ComputeElementDilatation(disp, elnum, disp_fespace, element_dilatation);
            (*dilatation)(elnum) = element_dilatation;
        }
    };

    void GlobalStressStrain::GlobalRotation(mfem::GridFunction *disp, mfem::GridFunction *rotation)
    {
        int numels = disp_fespace->GetNE();
        int dim = mesh->Dimension();
        // In 2D, rotation is a scalar. In 3D, it is a vector with 3 components.
        int rot_comp = (dim == 2) ? 1 : 3;

        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector element_rotation;
            ElementStressStrain Element;
            Element.ComputeElementRotation(disp, elnum, disp_fespace, element_rotation);
            for (int comp = 0; comp < rot_comp; comp++)
                (*rotation)(elnum + (numels * comp)) = element_rotation(comp);
        }
    };

    void GlobalStressStrain::GlobalStrainRotation(mfem::GridFunction *disp, mfem::GridFunction *strain, mfem::GridFunction *rotation)
    {
        int numels = L2_fespace->GetNE();
        int num_int_points = L2_fespace->GetFE(1)->GetNodes().GetNPoints();
        int dim = mesh->Dimension();
        // In 2D, there are 3 strain components. In 3D, there are 6 strain components.
        int str_comp = (dim == 2) ? 3 : 6;
        // In 2D, rotation is a scalar. In 3D, it is a vector with 3 components.
        int rot_comp = (dim == 2) ? 1 : 3;

        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector el_strain, el_rotation;
            ElementStressStrain Element;
            Element.ComputeElementStrainRotation(disp, elnum, disp_fespace, L2_fespace, el_strain, el_rotation);
            for (int ip = 0; ip < num_int_points; ip++)
            {
                for (int comp = 0; comp < str_comp; comp++)
                    (*strain)(num_int_points *elnum + ip + (num_int_points * numels * comp)) = el_strain(ip + comp * num_int_points);

                for (int comp = 0; comp < rot_comp; comp++)
                    (*rotation)(num_int_points *elnum + ip + (num_int_points * numels * comp)) = el_rotation(ip + comp * num_int_points);
            }
        }
    };

    void GlobalStressStrain::GlobalMaxShearStrainRotation(mfem::GridFunction *disp, mfem::GridFunction *max_strain, mfem::GridFunction *rotation)
    {
        int numels = L2_fespace->GetNE();
        int num_int_points = L2_fespace->GetFE(1)->GetNodes().GetNPoints();
        int dim = mesh->Dimension();
        // In 2D, rotation is a scalar. In 3D, it is a vector with 3 components.
        int rot_comp = (dim == 2) ? 1 : 3;

        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector el_strain, el_rotation;
            ElementStressStrain Element;
            Element.ComputeElementMaxShearStrainRotation(disp, elnum, disp_fespace, L2_fespace, el_strain, el_rotation);
            for (int ip = 0; ip < num_int_points; ip++)
            {
                (*max_strain)(num_int_points *elnum + ip) = el_strain(ip);

                for (int comp = 0; comp < rot_comp; comp++)
                    (*rotation)(num_int_points *elnum + ip + (num_int_points * numels * comp)) = el_rotation(ip + comp * num_int_points);
            }
        }
    }
    void GlobalStressStrain::GlobalStrainip(mfem::GridFunction &disp, mfem::GridFunction &strain)
    {
        int numels = L2_fespace->GetNE();
        int num_int_points = L2_fespace->GetFE(1)->GetNodes().GetNPoints();
        int dim = mesh->Dimension();
        // In 2D, there are 3 strain components. In 3D, there are 6 strain components.
        int str_comp = (dim == 2) ? 3 : 6;
        // In 2D, rotation is a scalar. In 3D, it is a vector with 3 components.
        int rot_comp = (dim == 2) ? 1 : 3;

        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector el_strain, el_rotation;
            ElementStressStrain Element;
            Element.ComputeElementStrainip(disp, elnum, disp_fespace, L2_fespace, el_strain);
            for (int ip = 0; ip < num_int_points; ip++)
            {
                for (int comp = 0; comp < str_comp; comp++)
                    strain(num_int_points * elnum + ip + (num_int_points * numels * comp)) = el_strain(ip + comp * num_int_points);
            }
        }
    }
}