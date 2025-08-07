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

#ifndef MFEM_CUSTOM_FUNCTIONS
#define MFEM_CUSTOM_FUNCTIONS

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

        // strain calculation with ParGridFunction
        void ComputeElementStrain(mfem::ParGridFunction &disp, int &elnum, mfem::ParFiniteElementSpace *parfes, mfem::Vector &elstrain);

        // stress calculation with E and NU
        void ComputeElementStress(mfem::Vector &elstrain, mfem::Coefficient &e, mfem::Coefficient &nu,
                                  int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress);

        // stress calculation with general stiffness matrix
        void ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat,
                                  int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &elstress);

        void ComputeElementStress(mfem::Vector &elstrain, mfem::Coefficient &e, mfem::Coefficient &nu,
                                  int &elnum, mfem::ParFiniteElementSpace *parfes, mfem::Vector &elstress);

        // // stress calculation with general stiffness matrix
        void ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat,
                                  int &elnum, mfem::ParFiniteElementSpace *parfes, mfem::Vector &elstress);

        void ComputeBoundaryElementArea(mfem::real_t &area, const mfem::FiniteElement *element, mfem::ElementTransformation *Trans);

        // dilatational and rotational strain and stress.
        void ComputeElementDilatation(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::real_t &dilatation);

        // strain calculation with ParGridFunction
        void ComputeElementRotation(mfem::GridFunction &disp, int &elnum, mfem::FiniteElementSpace *fes, mfem::Vector &rotation);
        ~ElementStressStrain() {};
    };

    class GlobalStressStrain
    {
    protected:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *fespace;
        mfem::ParMesh *pmesh;
        mfem::ParFiniteElementSpace *parfespace;

    private:
    public:
        GlobalStressStrain(mfem::Mesh *mesh_file, mfem::FiniteElementSpace *finite_el_space)
        {
            mesh = mesh_file;
            fespace = finite_el_space;
        };

        GlobalStressStrain(mfem::ParMesh *parmesh_file, mfem::ParFiniteElementSpace *parfinite_el_space)
        {
            pmesh = parmesh_file;
            parfespace = parfinite_el_space;
        };

        // Global strain calculation with GridFunction.
        void GlobalStrain(mfem::GridFunction &disp, mfem::GridFunction &strain);

        void GlobalStrain(mfem::ParGridFunction &disp, mfem::Array<mfem::ParGridFunction *> strain);

        // Global stress calculating with GridFunction.
        void GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &stress);

        void GlobalStress(mfem::GridFunction &strain, mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress);

        // Global stress calculating with ParGridFunction.
        void GlobalStress(mfem::Array<mfem::ParGridFunction *> strain, mfem::Coefficient &e,
                          mfem::Coefficient &nu, mfem::Array<mfem::ParGridFunction *> stress);

        void GlobalStress(mfem::Array<mfem::ParGridFunction *> strain, mfem::MatrixCoefficient &Cmat, mfem::Array<mfem::ParGridFunction *> stress);

        double ComputeBoundaryForce(mfem::GridFunction &stress, int &bdr_attribute, int &component);

        void GlobalDilatation(mfem::GridFunction &disp, mfem::GridFunction &dilatation);
        void GlobalRotation(mfem::GridFunction &disp, mfem::GridFunction &rotation);
        ~GlobalStressStrain() {};
    };
}

#endif