// This is a cpp file for custom operators written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
//------------------------------------------------------------------------------------------------------------
// 1. Construct operator to constrain displacement normal to boundary element node. The operator contains properly placed
// components of the orthonormal vector at each boundary node.
//------------------------------------------------------------------------------------------------------------

#include "customoperators.hpp"
#include "omp.h"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void ConstructNormalDisplacementConstraintOperator(mfem::FiniteElementSpace *fespace, mfem::Array<int> &constrained_boundary_elements, mfem::SparseMatrix *constraint_operator)
    {
        // Loop through the vector which contains the boundary element numbers of the constrained part of the boundary.
        // Extract the relevant boundary elements from the fespace, and from it the element transformation.
        // Set the integration points to be the computational nodes of the element and calculate the orthogonal vector at the node.
        // The orthogonal vector is NOT normalized.
        // In serial, the global dofs need to be determined for the local dofs of the element and then then the components of the normal
        // are to be placed in appropriate locations in the operator. The rest are zeros.
        // Sparse matrix height is the number of constrained nodes in the partition.
        // Spare matrix width is the total number of degrees of freedom in the partition.

        int dim = fespace->GetVDim();
        int num_constrained_els = constrained_boundary_elements.Size();

        for (int i = 0; i < num_constrained_els; i++)
        {
            int elnum = constrained_boundary_elements[i];
            const mfem::FiniteElement *element(fespace->GetBE(elnum));
            mfem::Array<int> eldofs;
            fespace->GetBdrElementDofs(elnum, eldofs);

            mfem::ElementTransformation *Trans(fespace->GetBdrElementTransformation(elnum));
            mfem::IntegrationRule int_rule(element->GetNodes());

            for (int node = 0; node < int_rule.Size(); node++)
            {
                mfem::Vector normal(dim), zero(dim);
                zero = 0.0;
                Trans->SetIntPoint(&(int_rule[node]));
                mfem::CalcOrtho(Trans->Jacobian(), normal);
                double normal_mag = normal.Norml2();
                add(zero, normal_mag, normal, normal);

                for (int vec_comp = 0; vec_comp < dim; vec_comp++)
                {
                    // Wrong row assignments. just constrained boundary dofs need to be assigned without repeats.
                    // Need to read through integration points.
                    constraint_operator->Set(eldofs[node], eldofs[node] + (vec_comp * fespace->GetNDofs()), normal(vec_comp));
                }
            };
        }
    };
}