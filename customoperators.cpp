// This is a cpp file for custom operators written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
//------------------------------------------------------------------------------------------------------------
// 1. Construct operator to constrain displacement normal to boundary element node. The operator contains properly placed
// components of the orthonormal vector at each boundary node.
//------------------------------------------------------------------------------------------------------------

#include "customoperators.hpp"
// #include "omp.h"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void RadialDisplacementConstraintOperator::Mult(const mfem::Vector &x, mfem::Vector &y) const
    {
        // This applies the (linear) left-hand side of the radial constraint for
        // every lateral-surface node:
        //     y_k = u1_k * X1_k + u2_k * X2_k
        y = 0.0;

        // Space is Ordering::byNODES, so component c of node j lives at
        // index j + c * num_nodes.

        MFEM_VERIFY(lateral_surface_dofs_x.Size() == y.Size(), "more surface nodes than constraint dofs");

        for (int k = 0; k < lateral_surface_dofs_x.Size(); k++)
        {
            y(k) = true_dof_coords(lateral_surface_dofs_x[k]) * x(lateral_surface_dofs_x[k]);
            y(k) += true_dof_coords(lateral_surface_dofs_y[k]) * x(lateral_surface_dofs_y[k]);
        }
    };

    // Transpose of Mult(): scatters the constraint multipliers back onto the
    // displacement dofs. If Mult maps u -> (C u), this maps lambda -> (C^T lambda):
    //
    //     y(j)             += X1_j * lambda_k
    //     y(j + num_nodes) += X2_j * lambda_k
    //
    // for each lateral-surface node j (with constraint index k). Used for the
    // (0,1) off-diagonal block of the saddle-point system.
    void RadialDisplacementConstraintOperator::MultTranspose(const mfem::Vector &x, mfem::Vector &y) const
    {
        y = 0.0;

        MFEM_VERIFY(lateral_surface_dofs_x.Size() == x.Size(), "more surface nodes than constraint dofs");

        for (int k = 0; k < lateral_surface_dofs_x.Size(); k++)
        {
            y(lateral_surface_dofs_x[k]) = true_dof_coords(lateral_surface_dofs_x[k]) * x(k);
            y(lateral_surface_dofs_y[k]) = true_dof_coords(lateral_surface_dofs_y[k]) * x(k);
        }
    };

    // Nonlinear right-hand side of the radial constraint, assembled as a plain
    // vector for the current displacement iterate x:
    //
    //     b_k = -(u1_k^2 + u2_k^2) / 2   for each lateral-surface node k
    //
    // The fixed-point iteration starts from x = 0 (so b = 0 on the first pass)
    // and calls this again with the updated iterate each sweep. The node loop
    // mirrors Mult() so the constraint-dof ordering in b matches y.
    void RadialDisplacementConstraintOperator::AssembleRHS(const mfem::Vector &x, mfem::Vector &b) const
    {
        b = 0.0;

        MFEM_VERIFY(lateral_surface_dofs_x.Size() == b.Size(), "more surface nodes than constraint dofs");

        for (int k = 0; k < lateral_surface_dofs_x.Size(); k++)
        {
            const double u1 = x(lateral_surface_dofs_x[k]);
            const double u2 = x(lateral_surface_dofs_y[k]);
            b(k) = -0.5 * ((u1 * u1) + (u2 * u2));
        }
    };
}