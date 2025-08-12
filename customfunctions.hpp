// This is a header file for custom functions written for the MFEM library.
// 1. Assign dirichlet boundary attributes.
// 2. Assign neumann boundary attributes.
// 3. Assign element attributes.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

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

}

#endif