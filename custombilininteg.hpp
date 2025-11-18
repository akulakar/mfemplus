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
    class IsotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {

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
        IsotropicElasticityIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu)
        {
            young_mod = &e;
            poisson_ratio = &nu;
        }
        IsotropicElasticityIntegrator() {};

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
    class AnisotropicElasticityIntegrator : public mfem::BilinearFormIntegrator
    {

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
        AnisotropicElasticityIntegrator(mfem::MatrixCoefficient &CMat)
        {
            stiffness = &CMat;
        }

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };

    /** Integrator for deviatoric component of isotropic linear elasticity using shear modulus mu.
        $$
          a(u,v) = (\mathrm{C}_{ijkl} u_{k,l} v_{i,j} )
        $$
        This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space. */
    class IsotropicElasticityDeviatoricIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *shear_mod;

    private:
#ifndef MFEM_THREAD_SAFE
        mfem::Vector shape, divshape;
        mfem::DenseMatrix dshape, gshape;
#endif

        // PA extension

        const mfem::DofToQuad *maps;        ///< Not owned
        const mfem::GeometricFactors *geom; ///< Not owned
        int vdim, ndofs;
        const mfem::FiniteElementSpace *fespace; ///< Not owned.

        std::unique_ptr<mfem::QuadratureSpace> q_space;
        /// Coefficients projected onto q_space
        std::unique_ptr<mfem::CoefficientVector> nu_quad;
        /// Workspace vector
        std::unique_ptr<mfem::QuadratureFunction> q_vec;

        /// Set up the quadrature space
        void SetUpQuadratureSpaceAndCoefficients(const mfem::FiniteElementSpace &fes);

    public:
        IsotropicElasticityDeviatoricIntegrator(mfem::Coefficient &mu)
        {
            shear_mod = &mu;
        }
        IsotropicElasticityDeviatoricIntegrator() {};

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    // Modified integrators from original MFEM.
    /** Class for local mass matrix assembling $a(u,v) := (Q u, v)$. However, using Gauss-Lobatto quadrature. */
    class GLMassIntegrator : public mfem::BilinearFormIntegrator
    {
    protected:
#ifndef MFEM_THREAD_SAFE
        mfem::Vector shape, te_shape;
#endif
        mfem::Coefficient *Q;
        // PA extension
        const mfem::FiniteElementSpace *fespace;
        mfem::Vector pa_data;
        const mfem::DofToQuad *maps;                 ///< Not owned
        const mfem::GeometricFactors *geom;          ///< Not owned
        const mfem::FaceGeometricFactors *face_geom; ///< Not owned
        int dim, ne, nq, dofs1D, quad1D;

        // void AssembleEA_(mfem::Vector &ea, const bool add);

    public:
        using ApplyKernelType = void (*)(const int, const mfem::Array<mfem::real_t> &,
                                         const mfem::Array<mfem::real_t> &, const mfem::Vector &,
                                         const mfem::Vector &, mfem::Vector &, const int, const int);

        using DiagonalKernelType = void (*)(const int, const mfem::Array<mfem::real_t> &,
                                            const mfem::Vector &, mfem::Vector &, const int,
                                            const int);

        MFEM_REGISTER_KERNELS(ApplyPAKernels, ApplyKernelType, (int, int, int));
        MFEM_REGISTER_KERNELS(DiagonalPAKernels, DiagonalKernelType, (int, int, int));
        struct Kernels
        {
            Kernels();
        };

    public:
        GLMassIntegrator(const mfem::IntegrationRule *ir = nullptr);

        /// Construct a mass integrator with coefficient q
        GLMassIntegrator(mfem::Coefficient &q, const mfem::IntegrationRule *ir = NULL);

        /** Given a particular Finite Element computes the element mass matrix
            elmat. */
        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Trans,
                                   mfem::DenseMatrix &elmat) override;

        // void AssembleMF(const mfem::FiniteElementSpace &fes) override;

        // using BilinearFormIntegrator::AssemblePA;
        // void AssemblePA(const mfem::FiniteElementSpace &fes) override;

        // void AssemblePABoundary(const mfem::FiniteElementSpace &fes) override;

        // void AssembleEA(const mfem::FiniteElementSpace &fes, mfem::Vector &emat,
        //                 const bool add) override;

        // virtual void AssembleEABoundary(const mfem::FiniteElementSpace &fes, mfem::Vector &emat,
        //                                 const bool add) override;

        // virtual void AssembleDiagonalPA(mfem::Vector &diag) override;

        // void AssembleDiagonalMF(mfem::Vector &diag) override;

        // void AddMultMF(const mfem::Vector &, mfem::Vector &) const override;

        // void AddMultPA(const mfem::Vector &, mfem::Vector &) const override;

        // void AddMultTransposePA(const mfem::Vector &, mfem::Vector &) const override;

        // static const mfem::IntegrationRule &GetRule(const mfem::FiniteElement &trial_fe,
        //                                             const mfem::FiniteElement &test_fe,
        //                                             const mfem::ElementTransformation &Trans);

        bool SupportsCeed() const override { return mfem::DeviceCanUseCeed(); }

        const mfem::Coefficient *GetCoefficient() const { return Q; }

        template <int DIM, int D1D, int Q1D>
        static void AddSpecialization()
        {
            ApplyPAKernels::Specialization<DIM, D1D, Q1D>::Add();
            DiagonalPAKernels::Specialization<DIM, D1D, Q1D>::Add();
        }

    protected:
        // const mfem::IntegrationRule *GetDefaultIntegrationRule(
        //     const mfem::FiniteElement &trial_fe,
        //     const mfem::FiniteElement &test_fe,
        //     const mfem::ElementTransformation &trans) const override
        // {
        //     return &GetRule(trial_fe, test_fe, trans);
        // }
    };

    /** Modified integrator for the bilinear form $a(u,v) := (Q u, v)$,
    where $u=(u_1,\dots,u_n)$ and $v=(v_1,\dots,v_n)$, $u_i$ and $v_i$ are defined
    by scalar FE through standard transformation. However, using Gauss-Lobatto quadrature rule */
    class GLVectorMassIntegrator : public mfem::BilinearFormIntegrator
    {
    private:
        int vdim;
        mfem::Vector shape, te_shape, vec;
        mfem::DenseMatrix partelmat;
        mfem::DenseMatrix mcoeff;
        int Q_order;

    protected:
        mfem::Coefficient *Q;
        mfem::VectorCoefficient *VQ;
        mfem::MatrixCoefficient *MQ;
        // PA extension
        mfem::Vector pa_data;
        const mfem::DofToQuad *maps;        ///< Not owned
        const mfem::GeometricFactors *geom; ///< Not owned
        int dim, ne, nq, dofs1D, quad1D;

    public:
        /// Construct an integrator with coefficient 1.0
        GLVectorMassIntegrator()
            : vdim(-1), Q_order(0), Q(NULL), VQ(NULL), MQ(NULL) {}
        /** Construct an integrator with scalar coefficient q.  If possible, save
            memory by using a scalar integrator since the resulting matrix is block
            diagonal with the same diagonal block repeated. */
        GLVectorMassIntegrator(mfem::Coefficient &q, int qo = 0)
            : vdim(-1), Q_order(qo), Q(&q), VQ(NULL), MQ(NULL) {}
        GLVectorMassIntegrator(mfem::Coefficient &q, const mfem::IntegrationRule *ir)
            : BilinearFormIntegrator(ir), vdim(-1), Q_order(0), Q(&q), VQ(NULL),
              MQ(NULL) {}
        /// Construct an integrator with diagonal coefficient q
        GLVectorMassIntegrator(mfem::VectorCoefficient &q, int qo = 0)
            : vdim(q.GetVDim()), Q_order(qo), Q(NULL), VQ(&q), MQ(NULL) {}
        /// Construct an integrator with matrix coefficient q
        GLVectorMassIntegrator(mfem::MatrixCoefficient &q, int qo = 0)
            : vdim(q.GetVDim()), Q_order(qo), Q(NULL), VQ(NULL), MQ(&q) {}

        int GetVDim() const { return vdim; }
        void SetVDim(int vdim_) { vdim = vdim_; }

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Trans,
                                   mfem::DenseMatrix &elmat) override;
    };
}
#endif