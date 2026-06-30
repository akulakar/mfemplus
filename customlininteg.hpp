// This is a header file for custom linear form integrators written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.

#ifndef MFEM_CUSTOM_LININTEG
#define MFEM_CUSTOM_LININTEG

#include "mfem.hpp"

namespace mfemplus
{
    // Integrators for variational fracture.
    /** New DomainLFIntegrator for damage force term.
        F(v) = \int_{\Omega} C_{ijkl} u_{k,l} u_{i,j} v
    **/
    class FractureDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Vector shape;
        mfem::Vector eldofdisp, eldofdamage;
        mfem::Array<int> eldofs;
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::Coefficient *volumetric_pressure;
        double viscosity_term;
        mfem::GridFunction *disp_gf;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        mfem::Vector CBu, Bu;
        mfem::Vector body_pressure;
        int planeApprox;
        mfem::DenseMatrix dshape, gshape;

    public:
        /// Constructs a domain integrator with a given Coefficient
        FractureDamageLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, double viscosity, mfem::GridFunction &disp, mfem::GridFunction &damage, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0, mfem::Coefficient *vol_pressure = nullptr)
            : young_mod(&e), poisson_ratio(&nu), viscosity_term(viscosity), disp_gf(&disp), damage_gf(&damage), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_pressure(vol_pressure) {};

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    class FractureHistoryVariableLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Vector shape;
        mfem::Vector eldofdisp, eldofdamage;
        mfem::Array<int> eldofs;
        int oa, ob;
        mfem::Coefficient *young_mod, *poisson_ratio, *viscosity_coeff;
        mfem::GridFunction *disp_gf, *dmg_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        mfem::Vector CBu, Bu;

    private:
        mfem::DenseMatrix dshape, gshape;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient
        FractureHistoryVariableLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp, mfem::FiniteElementSpace *disp_fespace, int a = 2, int b = 0)
            : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp), disp_fes(disp_fespace), oa(a), ob(b) {}

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    // Fracture specific integrators end here.

    // Eigenstrain body force integrator.
    // C_{ijkl} \epsilon^g_{kl}
    class EigenstrainBodyForceLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Vector shape;
        mfem::Vector eldofdisp, eldofdamage;
        mfem::Array<int> eldofs;
        int oa, ob;
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::VectorCoefficient *growth_strain_coeff;
        mfem::DenseMatrix C, B, BtC; // stiffness, strain-displacement
        mfem::Vector eps_g;          // eigenstrain
        int planeApprox;

    private:
        mfem::DenseMatrix dshape, gshape;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient
        EigenstrainBodyForceLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::VectorCoefficient &growthCoeff, int plane_approximation = 0)
            : young_mod(&e), poisson_ratio(&nu), growth_strain_coeff(&growthCoeff), planeApprox(plane_approximation) {}

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    class PressureBodyForceLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Coefficient *pressure;
        mfem::DenseMatrix B;      // stiffness, strain-displacement
        mfem::Vector body_stress; // volumetric pressure
        mfem::DenseMatrix dshape, gshape;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient
        PressureBodyForceLFIntegrator(mfem::Coefficient &pressure_coeff)
            : pressure(&pressure_coeff) {}

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    class PressureBodyForceDamageLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Coefficient *pressure;
        // mfem::Coefficient *Ey_coeff, *nu_coeff;
        mfem::DenseMatrix B;      // stiffness, strain-displacement
        mfem::Vector body_stress; // volumetric pressure
        mfem::DenseMatrix dshape, gshape;
        mfem::Vector shape;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage, elvec_input;
        mfem::real_t k_epsilon;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient

        PressureBodyForceDamageLFIntegrator(mfem::Coefficient &pressure_coeff, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace)
            : pressure(&pressure_coeff), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace) {};

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Linear form integrators for nonlinear elasticity begin here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    // Integrates the external tractions, body forces and internal nodal forces all in one loop.

    class StVenantKirchoffInternalForceLFIntegrator : public mfem::LinearFormIntegrator
    {
    protected:
        mfem::Coefficient *Ey_coeff, *nu_coeff;
        mfem::DenseMatrix C, BGradDisp, BNL; // stiffness, strain-displacement
        mfem::Vector Elg, F;                 // Green lagrange strain
        mfem::DenseMatrix dshape, gshape;
        mfem::Vector Gradu, Egl, S;
        mfem::GridFunction *disp_gf, *strain_gf, *stress_gf; // Green-Lagrange strain and 2nd Piola-Kirchoff stress.
        mfem::FiniteElementSpace *disp_fes, *L2_fes;         // H1 fes for displacement and L2 fes for strain and stress.
        mfem::Array<int> eldofs;                             // scalar for damage
        mfem::Vector eldofdisp, elvec_input, elstrain_ave, elstress_ave;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient

        StVenantKirchoffInternalForceLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::GridFunction &strain_gridfunc, mfem::GridFunction &stress_gridfunc, mfem::FiniteElementSpace *disp_fespace)
            : Ey_coeff(&e), nu_coeff(&nu), disp_gf(&disp_gridfunc), strain_gf(&strain_gridfunc), stress_gf(&stress_gridfunc), disp_fes(disp_fespace) {};
        // Fill in the average strain and stress for each element into the strain and stress gf;

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };
}
#endif