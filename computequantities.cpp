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
#include "omp.h"
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
    //------------------------------------------------------------------------------------------------------------------------------------

    void GlobalStressStrain::GlobalStrain(mfem::GridFunction &disp, mfem::GridFunction &strain)
    {
        int numels = fespace->GetNE();
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
            Element.ComputeElementStrain(disp, elnum, fespace, elstrain);
            for (int comp = 0; comp < str_comp; comp++)
                strain(elnum + (numels * comp)) = elstrain(comp);
        }
    };

    void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e,
                                          mfem::Coefficient &nu, mfem::GridFunction &stress)
    {

        int numels = fespace->GetNE();
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
            Element.ComputeElementStress(elstrain, e, nu, elnum, fespace, elstress);

            for (int comp = 0; comp < str_comp; comp++)
                stress(elnum + (numels * comp)) = elstress(comp);
        }
    };

    void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain,
                                          mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress)
    {

        int numels = fespace->GetNE();
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
            Element.ComputeElementStress(elstrain, Cmat, elnum, fespace, elstress);

            for (int comp = 0; comp < str_comp; comp++)
                stress(elnum + (numels * comp)) = elstress(comp);
        }
    };

    double GlobalStressStrain::ComputeBoundaryForce(mfem::GridFunction &stress, int &bdr_attribute, int &component)
    {
        double top_force = 0.0;

        int num_bdr_els = fespace->GetNBE();
        int num_els = fespace->GetNE();
        for (int i = 0; i < num_bdr_els; i++)
        {
            int bdr_attr = fespace->GetBdrAttribute(i);
            if (bdr_attr == bdr_attribute)
            {
                mfem::FaceElementTransformations *ftr = fespace->GetMesh()->GetBdrFaceTransformations(i);
                int volume_elem_index = ftr->Elem1No;
                double element_sig = stress(volume_elem_index + (component - 1) * num_els);
                double element_force = 0.0;
                mfem::real_t element_area;
                const mfem::FiniteElement *surf_el = fespace->GetBE(i);
                mfem::ElementTransformation *el_trans = fespace->GetBdrElementTransformation(i);
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
        int numels = fespace->GetNE();
        int dim = mesh->Dimension();
        // There is one measure of dilatation per element.
        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::real_t element_dilatation;
            ElementStressStrain Element;
            Element.ComputeElementDilatation(disp, elnum, fespace, element_dilatation);
            (*dilatation)(elnum) = element_dilatation;
        }
    };

    void GlobalStressStrain::GlobalRotation(mfem::GridFunction *disp, mfem::GridFunction *rotation)
    {
        int numels = fespace->GetNE();
        int dim = mesh->Dimension();
        // In 2D, distortion is a scalar. In 3D, it is a vector with 3 components.
        int dis_comp = (dim == 2) ? 1 : 3;

        // #pragma omp parallel for
        for (int elnum = 0; elnum < numels; elnum++)
        {
            mfem::Vector element_rotation;
            ElementStressStrain Element;
            Element.ComputeElementRotation(disp, elnum, fespace, element_rotation);
            for (int comp = 0; comp < dis_comp; comp++)
                (*rotation)(elnum + (numels * comp)) = element_rotation(comp);
        }
    };
}