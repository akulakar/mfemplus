// This is a header file for custom operators written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
//------------------------------------------------------------------------------------------------------------
// 1. Construct operator to constrain displacement normal to boundary element node. The operator contains properly placed
// components of the orthonormal vector at each boundary node.
//------------------------------------------------------------------------------------------------------------

#ifndef MFEM_CUSTOM_OPERATORS
#define MFEM_CUSTOM_OPERATORS

#include "mfem.hpp"
#include <memory>
#include <vector>

namespace mfemplus
{
    // This operator constrains the displacement field radially. The constraint on the displacement field for each point on the surface is:
    // u1 X1 + u2 X2 = -(u1^2 + u2^2)/2
    class RadialDisplacementConstraintOperator : public mfem::Operator
    {
    protected:
        mfem::GridFunction *node_coords;
        mfem::Vector true_dof_coords;
        mfem::Array<int> lateral_surface_dofs_x, lateral_surface_dofs_y;

    public:
        RadialDisplacementConstraintOperator(mfem::ParFiniteElementSpace *disp_fes, mfem::ParMesh *pmesh_in, mfem::Array<int> &lateral_surface_dofs_x_in, mfem::Array<int> &lateral_surface_dofs_y_in) : Operator(lateral_surface_dofs_x_in.Size(), disp_fes->GetTrueVSize()), lateral_surface_dofs_x(lateral_surface_dofs_x_in), lateral_surface_dofs_y(lateral_surface_dofs_y_in), node_coords(pmesh_in->GetNodes())
        {
            node_coords->GetTrueDofs(true_dof_coords);
        }

        void Mult(const mfem::Vector &x, mfem::Vector &y) const override;

        void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override;

        void AssembleRHS(const mfem::Vector &x, mfem::Vector &b) const;

        ~RadialDisplacementConstraintOperator() override
        {
        }
    };

}

#endif