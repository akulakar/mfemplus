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
        mfem::DenseMatrix C, B, CB, elmat_input;

        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

        // plane strain is 0, plane stress is 1.
        int planeApprox;

    public:
        IsotropicElasticityIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), planeApprox(plane_approximation)
        {
            // int dof = el.GetDof();
            // int dim = el.GetDim();
            // int str_comp = (dim == 2) ? 3 : 6;
            // C.SetSize(str_comp, str_comp);  // Stiffness in Voigt form
            // B.SetSize(str_comp, dof * dim); // Strain displacement matrix
            // CB.SetSize(str_comp, dof * dim);
            // elmat_input.SetSize(dof * dim, dof * dim);
        }

        void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::DenseMatrix &elmat) override;
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
        mfem::DenseMatrix C, B, CB, elmat_input;

        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        AnisotropicElasticityIntegrator(mfem::MatrixCoefficient &CMat) : stiffness(&CMat) {}

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

        mfem::Vector shape, divshape;
        mfem::DenseMatrix dshape, gshape;
        mfem::DenseMatrix B, B_vol, B_dev, elmat_input; // Strain displacement matrix
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
        GLMassIntegrator(const mfem::IntegrationRule *ir)
            : mfem::BilinearFormIntegrator(ir), Q(nullptr), maps(nullptr), geom(nullptr)
        {
            // static GLMassIntegrator::Kernels kernels;
        }

        // Construct a mass integrator with coefficient q
        GLMassIntegrator(mfem::Coefficient &q, const mfem::IntegrationRule *ir)
            : GLMassIntegrator(ir)
        {
            Q = &q;
        }

        /** Given a particular Finite Element computes the element mass matrix
            elmat. */
        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Trans,
                                   mfem::DenseMatrix &elmat) override;

        bool SupportsCeed() const override { return mfem::DeviceCanUseCeed(); }

        const mfem::Coefficient *GetCoefficient() const { return Q; }

        template <int DIM, int D1D, int Q1D>
        static void AddSpecialization()
        {
            ApplyPAKernels::Specialization<DIM, D1D, Q1D>::Add();
            DiagonalPAKernels::Specialization<DIM, D1D, Q1D>::Add();
        }
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
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Integrators for variational fracture.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    /** Integrator for variational fracture strain energy term using general stiffness tensor, specifying Young's modulus E and Poisson's ratio nu, and multiplied with the degradation damage function (1 - d^{2}) + k_{\epsilon}.
       $$
         a(u,v) =  ((1 - d)^{2} + k_{\epsilon})*(\mathrm{C}_{ijkl} u_{k,l} v_{i,j})
       $$
       This is a 'Vector' integrator, i.e. defined for FE spaces using multiple copies of a scalar FE space. **/
    class IsotropicElasticityDamageIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::real_t k_epsilon;
        mfem::GridFunction *damage_gf;
        mfem::FiniteElementSpace *damage_fes;
        mfem::Array<int> eldofs; // scalar for damage
        mfem::Vector eldofdamage;
        mfem::DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        int planeApprox;
        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        IsotropicElasticityDamageIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::real_t k_eps, mfem::GridFunction &damage_gridfunc, mfem::FiniteElementSpace *damage_fespace, int plane_approximation = 0) : young_mod(&e), poisson_ratio(&nu), k_epsilon(k_eps), damage_gf(&damage_gridfunc), damage_fes(damage_fespace), planeApprox(plane_approximation) {};

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };

    /** Integrator for variational fracture, elastic energy and damage term, using Young's modulus E and Poisson's ratio nu
       $$
         a(u,v) =  (\mathrm{C}_{ijkl} u_{k,l} u_{i,j}) d \phi
       $$
       This is a scalar integrator, i.e. defined for scalar FE spaces. **/
    class IsotropicStrainEnergyDamageIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::Coefficient *volumetric_pressure;
        mfem::GridFunction *disp_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::Array<int> eldofs;
        mfem::Vector eldofdamage;
        mfem::Vector eldofdisp;
        mfem::Vector CBu, Bu;
        mfem::DenseMatrix C, B, CB, elmat_input; // Stiffness in Voigt form
        mfem::Vector body_pressure;

        int planeApprox, GershgorinCheck;

        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        IsotropicStrainEnergyDamageIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::FiniteElementSpace *disp_fespace, int plane_approximation = 0, mfem::Coefficient *vol_press = nullptr) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace), planeApprox(plane_approximation), volumetric_pressure(vol_press) {};

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Fracture integrators end here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    // Nonlinear elasticity integrators declared here.
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    /**
        Integrator for St. Venant Kirchoff tangent stiffness matrix
       $$
         a(\Delta u,v) =  \int_{\Omega_{e}} \bar{E}(\Delta u, v) : C : \Delta E(u, \Delta u) + C : E(u) : \Delta \bar{E} (\Delta u, v) ~\mathrm{d}\Omega_{e}
       $$
       This is a vector integrator, i.e. defined for vector FE spaces.
    **/
    class StVenantKirchoffTangentStiffnessIntegrator : public mfem::BilinearFormIntegrator
    {

    protected:
        mfem::Coefficient *young_mod, *poisson_ratio;
        mfem::GridFunction *disp_gf;
        mfem::FiniteElementSpace *disp_fes;
        mfem::Array<int> eldofs;
        mfem::Vector eldofdisp;
        mfem::Vector CBu, Bu, F, Egl, Gradu, S;
        mfem::DenseMatrix C, BGradDisp, BNL, CB, Sigma, elmat_input_1, elmat_input_2; // Stiffness in Voigt form, Bmatrix (for displacement gradients);

        mfem::Vector shape;
        mfem::DenseMatrix dshape, gshape;

    public:
        StVenantKirchoffTangentStiffnessIntegrator(mfem::Coefficient &e, mfem::Coefficient &nu, mfem::GridFunction &disp_gridfunc, mfem::ParFiniteElementSpace *disp_fespace) : young_mod(&e), poisson_ratio(&nu), disp_gf(&disp_gridfunc), disp_fes(disp_fespace) {};

        void AssembleElementMatrix(const mfem::FiniteElement &el,
                                   mfem::ElementTransformation &Tr,
                                   mfem::DenseMatrix &elmat) override;
    };
}
#endif