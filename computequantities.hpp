// This is a header file for custom functions written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
// List of functions are given below.
// 1. Elasticity
//      (a) Recover strain vector for each element using B matrix given a displacement grid function.
//      (b) Recover stress vector for each element given strain and constitutive law.
//      (c) Recover dilatational strain and stress.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// 2.

#ifndef MFEM_COMPUTE_QUANTITIES
#define MFEM_COMPUTE_QUANTITIES

#include "mfem.hpp"
#include <memory>
#include <vector>

namespace mfemplus
{
    class AccessMFEMFunctions : public mfem::BilinearFormIntegrator
    {
    public:
        using mfem::BilinearFormIntegrator::GetIntegrationRule;
    };

    struct ElementStressStrain
    {
        friend class GlobalStressStrain;

    public:
    protected:
    private:
        // strain calculation with GridFunction.
        void ComputeElementStrain(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstrain);

        // stress calculation with E and NU
        void ComputeElementStress(mfem::Vector &elstrain, mfem::Coefficient &e, mfem::Coefficient &nu,
                                  int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress);

        // stress calculation with general stiffness matrix
        void ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat,
                                  int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress);

        void ComputeBoundaryElementArea(mfem::real_t &area, const mfem::FiniteElement *element, mfem::ElementTransformation *Trans);

        // dilatational and rotational strain and stress.
        void ComputeElementDilatation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::real_t &dilatation);

        void ComputeElementRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &rotation);
        void ComputeElementRotationip(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &rotation);
        void ComputeElementMaxShearStrainip(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &shearstrain);

        mfem::Vector ComputeElementDisplacementGradients(mfem::DenseMatrix &gshape, mfem::Vector &eldofdisp);
        void ComputeElementStrainip(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &elstrain);
        void ComputeElementStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &strain, mfem::Vector &rotation);
        void ComputeElementMaxShearStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &max_shear_strain, mfem::Vector &rotation);
        void ComputeElementStrainMaxShearStrainRotation(mfem::GridFunction *disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &strain, mfem::Vector &max_shear_strain, mfem::Vector &rotation);
        void ComputeElementAverageMaxShearPrincipalStrainRotation(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *disp_fes, mfem::FiniteElementSpace *L2_fes, mfem::Vector &principal_strains, double &max_shear_strain, mfem::Vector &rotation);

        ~ElementStressStrain() {};
    };

    class GlobalStressStrain
    {
    protected:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *disp_fespace, *L2_fespace;
        // mfem::ParMesh *pmesh;
        // mfem::ParFiniteElementSpace *parfespace;

    private:
    public:
        GlobalStressStrain(mfem::Mesh *mesh_file, mfem::FiniteElementSpace *disp_finite_el_space, mfem::FiniteElementSpace *L2_finite_el_space)
        {
            mesh = mesh_file;
            disp_fespace = disp_finite_el_space;
            L2_fespace = L2_finite_el_space;
        };
        // Global strain calculation with GridFunction.
        void GlobalStrain(mfem::GridFunction &disp, mfem::GridFunction &strain);

        // Global stress calculating with GridFunction.
        void GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &stress);

        void GlobalStress(mfem::GridFunction &strain, mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress);

        double ComputeBoundaryForce(mfem::GridFunction &stress, int &bdr_attribute, int &component);

        // Global dilatation and rotation
        void GlobalDilatation(mfem::GridFunction *disp, mfem::GridFunction *dilatation);
        void GlobalRotation(mfem::GridFunction *disp, mfem::GridFunction *rotation);
        void GlobalRotationip(mfem::GridFunction &disp, mfem::GridFunction &rotation);
        void GlobalStrainip(mfem::GridFunction &disp, mfem::GridFunction &strain);
        void GlobalMaxShearStrainip(mfem::GridFunction &disp, mfem::GridFunction &maxshearstrain);

        void GlobalStrainRotation(mfem::GridFunction *disp, mfem::GridFunction *strain, mfem::GridFunction *rotation);
        void GlobalMaxShearStrainRotation(mfem::GridFunction *disp, mfem::GridFunction *max_strain, mfem::GridFunction *rotation);
        void GlobalAverageMaxShearPrincipalStrainRotation(mfem::GridFunction &disp, mfem::GridFunction &principal_strains, mfem::GridFunction &max_shear_strain, mfem::GridFunction &rotation);
        ~GlobalStressStrain() {};
    };
}

#endif