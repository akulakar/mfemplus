//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// This is a header file for custom implementations of hyperbolic conservation laws written for the MFEM library.
// This file contains general hyperbolic conservation element/face form
// integrators. HyperbolicFormIntegrator and NumericalFlux are defined.
//
// HyperbolicFormIntegrator is a NonlinearFormIntegrator that implements
// element weak divergence and interface flux
//
//     ∫_K F(u):∇v,   -∫_f F̂(u)⋅n[v]
//
// Here, K is an element, f is a face, n normal and [⋅] is jump. This form
// integrator is coupled with NumericalFlux that implements the numerical flux
// F̂. The following numerical fluxes have been implemented
//
// To implement a specific hyperbolic conservation laws, users can create
// derived classes from FluxFunction with overloaded ComputeFlux. One can
// optionally overload ComputeFluxDotN to avoid creating dense matrix when
// computing normal flux.
//
// At each call of HyperbolicFormIntegrator::AssembleElementVector
// HyperbolicFormIntegrator::AssembleFaceVector, the maximum characteristic
// speed will be updated. This will not be reinitialized automatically. To
// reinitialize, use HyperbolicFormIntegrator::ResetMaxCharSpeed.
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#ifndef MFEM_CUSTOM_HYPERBOLIC
#define MFEM_CUSTOM_HYPERBOLIC

#include "mfem.hpp"
#include <memory>
#include <vector>

namespace mfemplus
{
    // @brief
    // The nonlinearForm will only evaluate boundary face integrators at the marked boundaries.
    // However, we need to provide the state vector for the inflow condition at each marked boundary.
    class HyperbolicFormBdrInflowIntegrator : public mfem::NonlinearFormIntegrator
    {
    protected:
        // The maximum characteristic speed, updated during element/face vector assembly
        mfem::real_t max_char_speed;
        const mfem::NumericalFlux &numFlux; // Numerical flux that maps F(u±,x) to F̂
        const mfem::FluxFunction &fluxFunction;
        const int IntOrderOffset; // integration order offset, 2*p + IntOrderOffset.
        const mfem::real_t sign;

        // Local storage for element integration

        mfem::Vector GhostBdrState; // Prescribed boundary conditions will act as a ghost state for a flux
        mfem::Vector ElementState;  // state value at an integration point in the element in the mesh
        mfem::Vector shape;         // shape function value at an integration point for the element in the mesh
        mfem::DenseMatrix flux;     // flux value at an integration point

        mfem::Vector nor;   // normal vector, see mfem::CalcOrtho()
        mfem::Vector fluxN; // F̂(u±,x) n

    public:
        const int num_equations; // the number of equations
        /**
         * @brief Construct a new Hyperbolic Form Integrator object for inflow at the boundary with ghost states.
         *
         * @param[in] numFlux numerical flux
         * @param[in] IntOrderOffset integration order offset
         * @param[in] sign sign of the convection term
         */
        HyperbolicFormBdrInflowIntegrator(const mfem::NumericalFlux &numFlux, mfem::Vector &GhostBdrState_in, const int IntOrderOffset = 0, const mfem::real_t sign = 1.)
            : NonlinearFormIntegrator(),
              numFlux(numFlux),
              fluxFunction(numFlux.GetFluxFunction()),
              GhostBdrState(GhostBdrState_in),
              IntOrderOffset(IntOrderOffset), sign(sign),
              num_equations(fluxFunction.num_equations)
        {
            flux.SetSize(num_equations, fluxFunction.dim);
            GhostBdrState.SetSize(num_equations);
            ElementState.SetSize(num_equations);
            fluxN.SetSize(num_equations);
            nor.SetSize(fluxFunction.dim);
            ResetMaxCharSpeed();
        };

        /**
         * @brief Reset the Max Char Speed 0
         *
         */
        void ResetMaxCharSpeed()
        {
            max_char_speed = 0.0;
        }

        mfem::real_t GetMaxCharSpeed()
        {
            return max_char_speed;
        }

        /**
         * @brief Implements <-F̂(u⁻,u⁺,x) n, [v]> with abstract F̂ computed by
         * NumericalFlux::Eval() of the numerical flux object
         *
         * @param[in] el1 finite element of the first element
         * @param[in] el2 finite element of the second element
         * @param[in] Tr face element transformations
         * @param[in] elfun local coefficient of basis from both elements
         * @param[out] elvect evaluated dual vector <-F̂(u⁻,u⁺,x) n, [v]>
         */
        void AssembleFaceVector(const mfem::FiniteElement &el1,
                                const mfem::FiniteElement &el2,
                                mfem::FaceElementTransformations &Tr,
                                const mfem::Vector &elfun, mfem::Vector &elvect) override;
    };
}

#endif