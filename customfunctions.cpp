// This is a cpp file for custom functions written for the MFEM library.
// 1. Assign dirichlet boundary attributes.
// 2. Assign neumann boundary attributes.
// 3. Assign element attributes.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
#include "customfunctions.hpp"
// #include "omp.h"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void AssignAttributes(mfem::Mesh *mesh)
    {
        int num_els = mesh->GetNE();
        int num_bdr_els = mesh->GetNBE();

        for (int elnum = 0; elnum < num_els; elnum++)
            mesh->SetAttribute(elnum, elnum + 1);

        for (int bdr_elnum = 0; bdr_elnum < num_bdr_els; bdr_elnum++)
            mesh->SetBdrAttribute(bdr_elnum, bdr_elnum + 1);

        mesh->SetAttributes();
    };
}