// This is a header file for custom integrators written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
// List of integrators is given below.
// 1. Elasticity
//      (a) Isotropic elasticity with E and nu.
//      (b) Fully anisotropic elasticty with 21 independent constants.
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// 2.

#ifndef MFEM_CUSTOM_BILININTEG
#define MFEM_CUSTOM_BILININTEG

#include "mfem.hpp"

namespace mfemplus
{
    /** Integrator for linear elasticity using general stiffness tensor, but specifying Young's modulus E and Poisson's ratio nu:
        $$
          a(u,v) = (\mathrm{C}_{ijkl} u_{k,l} v_{i,j} )
        $$
        This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space. */
    class ThreeDIsotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {
        friend class ElasticityComponentIntegrator;

    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;

    private:
#ifndef MFEM_THREAD_SAFE
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;
#endif

        // PA extension

        const mfem::DofToQuad *maps;        ///< Not owned
        const mfem::GeometricFactors *geom; ///< Not owned
        int vdim, ndofs;
        const mfem::FiniteElementSpace *fespace; ///< Not owned.

        std::unique_ptr<mfem::QuadratureSpace> q_space;
        /// Coefficients projected onto q_space
        std::unique_ptr<mfem::CoefficientVector> E_quad, nu_quad;
        /// Workspace vector
        std::unique_ptr<mfem::QuadratureFunction> q_vec;

        /// Set up the quadrature space
        void SetUpQuadratureSpaceAndCoefficients(const mfem::FiniteElementSpace &fes);

    public:
        ThreeDIsotropicElasticityIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu)
        {
            young_mod = &e;
            poisson_ratio = &nu;
        }
        ThreeDIsotropicElasticityIntegrator() {};

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
        };

    /** Integrator for anistotropic linear elasticity using general stiffness tensor. A total of 21 elastic constants need to be prescribed:
        $$
          a(u,v) = ( \mathrm{C}_{ijkl} u_{k,l} v_{i,j} )
          a(u,v) = ( \boldsymbol{\epsilon(v)}^T \boldsymbol{B} \boldsymbol{\epsilon(u)} )
        $$
        This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space. */
    class ThreeDAnisotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {
        friend class ElasticityComponentIntegrator;

    protected:
        mfem::MatrixCoefficient *stiffness;

    private:
#ifndef MFEM_THREAD_SAFE
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;
#endif

        // PA extension

        const mfem::DofToQuad *maps;        ///< Not owned
        const mfem::GeometricFactors *geom; ///< Not owned
        int vdim, ndofs;
        const mfem::FiniteElementSpace *fespace; ///< Not owned.

        std::unique_ptr<mfem::QuadratureSpace> q_space;
        /// Workspace vector
        std::unique_ptr<mfem::QuadratureFunction> q_vec;

        /// Set up the quadrature space
        void SetUpQuadratureSpaceAndCoefficients(const mfem::FiniteElementSpace &fes);

    public:
        ThreeDAnisotropicElasticityIntegrator(mfem::MatrixCoefficient &CMat)
        {
            stiffness = &CMat;
        }

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };
}

#endif