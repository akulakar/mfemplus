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

    void ConstructNormalDisplacementConstraintOperator(mfem::FiniteElementSpace *fespace, mfem::Array<int> &constrained_boundary_elements, mfem::SparseMatrix *constraint_operator);
}

#endif