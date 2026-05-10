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
        int oa, ob;
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::Coefficient *viscosity_term;
        mfem::GridFunction *disp_gf;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::DenseMatrix C, B, CB; // stiffness, strain-displacement, Stiffness times strain-displacement in Voigt form
        mfem::Vector CBu, Bu;

    private:
        mfem::DenseMatrix dshape, gshape;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient
        FractureDamageLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::Coefficient &viscosity, mfem::GridFunction &disp, mfem::GridFunction &damage, mfem::FiniteElementSpace *disp_fespace, int a = 2, int b = 0)
            : young_mod(&e), poisson_ratio(&nu), viscosity_term(&viscosity), disp_gf(&disp), damage_gf(&damage), disp_fes(disp_fespace), oa(a), ob(b) {}

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

    private:
        mfem::DenseMatrix dshape, gshape;

    public:
        // Need to modify the constructors...

        /// Constructs a domain integrator with a given Coefficient
        EigenstrainBodyForceLFIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::VectorCoefficient &growthCoeff, int a = 2, int b = 0)
            : young_mod(&e), poisson_ratio(&nu), growth_strain_coeff(&growthCoeff), oa(a), ob(b) {}

        void AssembleDevice(const mfem::FiniteElementSpace &fes, const mfem::Array<int> &markers, mfem::Vector &b) override {};

        /** Given a particular Finite Element and a transformation (Tr)
            computes the element right hand side element vector, elvect. **/
        void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect) override;

        virtual void AssembleRHSElementVect(const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr, mfem::Vector &elvect) override {};

        using mfem::LinearFormIntegrator::AssembleRHSElementVect;
    };
}
#endif