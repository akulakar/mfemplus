// This is a header file for custom functions written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
// List of functions are given below.
// 1. Elasticity
//      (a) Recover strain tensor for each element using B matrix given a displacement grid function.
//      (b) Recover stress tensor for each element given strain and constitutive law.
//
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

    class ElementStressStrain : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *elastic_mod, *poisson_ratio;
        mfem::MatrixCoefficient *elastic_constants;
        const mfem::FiniteElement *el;
        mfem::ElementTransformation *Tr;
        mfem::FiniteElementSpace *fespace;
        mfem::ParFiniteElementSpace *parfespace;
        const int *elnumber;

    private:
#ifndef MFEM_THREAD_SAFE
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;
#endif
    public:
        ElementStressStrain(const mfem::FiniteElement &element, mfem::ElementTransformation &Trans,
                            const int &elnum, mfem::FiniteElementSpace *fes)
        {
            el = &element;
            Tr = &Trans;
            elnumber = &elnum;
            fespace = fes;
        }

        ElementStressStrain(const mfem::FiniteElement &element, mfem::ElementTransformation &Trans,
                            const int &elnum, mfem::ParFiniteElementSpace *parfes)
        {
            el = &element;
            Tr = &Trans;
            elnumber = &elnum;
            parfespace = parfes;
        }

        ElementStressStrain(const mfem::FiniteElement &element, mfem::ElementTransformation &Trans)
        {
            el = &element;
            Tr = &Trans;
        }

        // strain calculation with GridFunction.
        void ComputeElementStrain(mfem::GridFunction &disp, mfem::Vector &elstrain);

        // strain calculation with ParGridFunction
        void ComputeElementStrain(mfem::ParGridFunction &disp, mfem::Vector &elstrain);

        void ComputeElementStress(mfem::Vector &elstrain, mfem::Coefficient &e,
                                  mfem::Coefficient &nu, mfem::Vector &elstress);

        void ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat,
                                  mfem::Vector &elstress);

        void ComputeBoundaryElementArea(mfem::real_t &area);
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

        // Global strain calculation with ParGridFunction.
        void GlobalStrain(mfem::ParGridFunction &disp, mfem::ParGridFunction &strain);

        // Global stress calculating with GridFunction.
        void GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &stress);

        void GlobalStress(mfem::GridFunction &strain, mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress);

        // Global stress calculating with ParGridFunction.
        void GlobalStress(mfem::ParGridFunction &strain, mfem::Coefficient &e, mfem::Coefficient &nu, mfem::ParGridFunction &stress);

        void GlobalStress(mfem::ParGridFunction &strain, mfem::MatrixCoefficient &Cmat, mfem::ParGridFunction &stress);


        void Global3DStrainComponents(mfem::GridFunction &strain, mfem::GridFunction &eps11, mfem::GridFunction &eps22,
                                      mfem::GridFunction &eps33, mfem::GridFunction &eps23, mfem::GridFunction &eps13,
                                      mfem::GridFunction &eps12);

        void Global3DStressComponents(mfem::GridFunction &stress, mfem::GridFunction &sig11, mfem::GridFunction &sig22,
                                      mfem::GridFunction &sig33, mfem::GridFunction &sig23, mfem::GridFunction &sig13,
                                      mfem::GridFunction &sig12);

        // std::vector<std::unique_ptr<mfem::GridFunction>> straincomponents(mfem::GridFunction &strain);
        // std::vector<std::unique_ptr<mfem::GridFunction>> stresscomponents(mfem::GridFunction &stress);
    };
}

#endif