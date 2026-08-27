//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// This is a cpp file for custom implementations of hyperbolic conservation laws written for the MFEM library.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
#include "customhyperbolic.hpp"
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
namespace mfemplus
{
    void HyperbolicFormBdrInflowIntegrator::AssembleFaceVector(
        const mfem::FiniteElement &el1, const mfem::FiniteElement &el2,
        mfem::FaceElementTransformations &Tr, const mfem::Vector &elfun, mfem::Vector &elvect)
    {
        // This integrator is only for boundary faces, so both el1 and el2 are the same! We will deal with just el1.
        // First, check the boundary attribute, and if it is marked for inflow, then evaluate the corresponding state.

        // current elements' the number of degrees of freedom
        // does not consider the number of equations
        const int dof1 = el1.GetDof();

        shape.SetSize(dof1);

        elvect.SetSize(dof1 * num_equations);
        elvect = 0.0;

        const mfem::DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equations);
        mfem::DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equations);

        // Obtain integration rule. If integration is rule is given, then use it.
        // Otherwise, get (2*p + IntOrderOffset) order integration rule
        const mfem::IntegrationRule *ir = IntRule;
        if (!ir)
        {
            const int order = 2 * std::max(el1.GetOrder(), el2.GetOrder()) + IntOrderOffset;
            ir = &mfem::IntRules.Get(Tr.GetGeometryType(), order);
        }
        // loop over integration points
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            Tr.SetAllIntPoints(&ip); // set face and element int. points

            // Calculate basis functions on internal elements at the face
            el1.CalcShape(Tr.GetElement1IntPoint(), shape);

            // Interpolate elfun at the integration point
            elfun1_mat.MultTranspose(shape, ElementState);

            // Get the normal vector and the flux on the face
            if (nor.Size() == 1) // if 1D, use 1 or -1.
            {
                // This assume the 1D integration point is in (0,1). This may not work
                // if this changes.
                nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;
            }
            else
            {
                CalcOrtho(Tr.Jacobian(), nor);
            }
            // Compute F(u+, x) and F(u-, x) with maximum characteristic speed
            // Compute hat(F) using evaluated quantities
            const mfem::real_t speed = numFlux.Eval(GhostBdrState, ElementState, nor, Tr, fluxN);

            // Update the global max char speed
            max_char_speed = std::max(speed, max_char_speed);

            // pre-multiply integration weight to flux
            AddMult_a_VWt(ip.weight * sign, shape, fluxN, elvect1_mat);
            // AddMult_a_VWt(+ip.weight * sign, shape2, fluxN, elvect2_mat);
        }
    }
}