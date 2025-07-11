// This constains the implementation of custom functions written for the MFEM library.
// Written by members of the Applied Mechanics Lab at Brown university.
// List of functions are given below.
// 1. Elasticity
//      (a) Recover strain tensor for each element using B matrix given a displacement grid function.
//      (b) Recover stress tensor for each element given strain and constitutive law.
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


    #include "customfunctions.hpp"
    #include <cmath>
    #include <algorithm>
    #include <memory>

    using namespace std;
    namespace mfemplus
    {
        void ElementStressStrain::ComputeElementStrain(mfem::GridFunction &disp, mfem::Vector &elstrain)
        {

            // To recover element strain, we need the B matrix at each integration point multiplied by
            // the displacement vector at each node in the element. Therefore, the strain varies
            // within the element if the B matrix is not constant. However, the strain is averaged
            // within an element.

            int dof = el->GetDof();
            int dim = el->GetDim();

            mfem::Array<int> eldofs(dof * dim);
            mfem::Array<double> eldofdisp(dof * dim);

            fespace->GetElementVDofs(*elnumber, eldofs);
            int eltype = fespace->GetElementType(*elnumber);

            for (int i = 0; i < eldofdisp.Size(); i++){
                int dof = eldofs[i];
                eldofdisp[i] = disp(dof);
            }

            mfem::real_t w, E, NU;

            MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

    #ifdef MFEM_THREAD_SAFE
            DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
            Vector divshape(dim * dof);
    #else
            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);

    #endif
            if (dim == 2)
                elstrain.SetSize(3);
            else if (dim == 3)
                elstrain.SetSize(6);
            
            elstrain = 0.0;

            const mfem::IntegrationRule *ir = GetIntegrationRule(*el, *Tr);
            if (ir == NULL)
            {
                int order = 2 * Tr->OrderGrad(el); // correct order?
                ir = &mfem::IntRules.Get(el->GetGeomType(), order);
            }

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                el->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.

                Tr->SetIntPoint(&ip);
                mfem::Mult(dshape, Tr->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

                // Here we want to use voigt notation to compute the strain vector. 
                // The strain vector has 3 components in 2D and 6 components in 3D.
                mfem::DenseMatrix B;

                if (dim == 2)
                {
                    B.SetSize(3, dof * dim); // In 2D, we have 3 unique strain components.
                    B = 0.0;

                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf) = gshape(spf, 1);
                        B(2, spf + dof) = gshape(spf, 0);
                
                    }
                }

                else if (dim == 3)
                {
                    B.SetSize(6, dof * dim); // In 3D, we have 6 unique strain components.
                    B = 0.0;

                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf + 2 * dof) = gshape(spf, 2);
                        B(3, spf + dof) = gshape(spf, 2);
                        B(3, spf + 2 * dof) = gshape(spf, 1);
                        B(4, spf) = gshape(spf, 2);
                        B(4, spf + 2 * dof) = gshape(spf, 0);
                        B(5, spf) = gshape(spf, 1);
                        B(5, spf + dof) = gshape(spf, 0);
                    }
                }

            mfem::Vector temp(elstrain.Size()); temp = 0.0;
            B.Mult(eldofdisp,temp);
            double w = 1.0/(ir->GetNPoints());
            elstrain.Add(w,temp);
        }
    };

    void ElementStressStrain::ComputeElementStrain(mfem::ParGridFunction &disp, mfem::Vector &elstrain)
        { 

            // To recover element strain, we need the B matrix at each integration point multiplied by
            // the displacement vector at each node in the element. Therefore, the strain varies
            // within the element if the B matrix is not constant. However, the strain is averaged
            // within an element.

            int dof = el->GetDof();
            int dim = el->GetDim();

            mfem::Array<int> eldofs(dof * dim);
            mfem::Array<double> eldofdisp(dof * dim);

            parfespace->GetElementVDofs(*elnumber, eldofs);
            int eltype = parfespace->GetElementType(*elnumber);

            for (int i = 0; i < eldofdisp.Size(); i++){
                int dof = eldofs[i];
                eldofdisp[i] = disp(dof);
            }

            mfem::real_t w, E, NU;

            MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

    #ifdef MFEM_THREAD_SAFE
            DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
            Vector divshape(dim * dof);
    #else
            dshape.SetSize(dof, dim);
            gshape.SetSize(dof, dim);

    #endif
            if (dim == 2)
                elstrain.SetSize(3);
            else if (dim == 3)
                elstrain.SetSize(6);
            elstrain = 0.0;

            const mfem::IntegrationRule *ir = GetIntegrationRule(*el, *Tr);
        if (ir == NULL)
        {
            int order = 2 * Tr->OrderGrad(el); // correct order?
            ir = &mfem::IntRules.Get(el->GetGeomType(), order);
        }

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const mfem::IntegrationPoint &ip = ir->IntPoint(i);

            el->CalcDShape(ip, dshape); // Gradients of the shape functions in the reference element.

            Tr->SetIntPoint(&ip);
            mfem::Mult(dshape, Tr->InverseJacobian(), gshape); // Recovering the gradients of the shape functions in the physical space.

           // Here we want to use voigt notation to compute the strain vector. 
                // The strain vector has 3 components in 2D and 6 components in 3D.
                mfem::DenseMatrix B;

                if (dim == 2)
                {
                    B.SetSize(3, dof * dim); // In 2D, we have 3 unique strain components.
                    B = 0.0;

                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf) = gshape(spf, 1);
                        B(2, spf + dof) = gshape(spf, 0);
                
                    }
                }

                else if (dim == 3)
                {
                    B.SetSize(6, dof * dim); // In 3D, we have 6 unique strain components.
                    B = 0.0;

                    for (int spf = 0; spf < dof; spf++)
                    {
                        B(0, spf) = gshape(spf, 0);
                        B(1, spf + dof) = gshape(spf, 1);
                        B(2, spf + 2 * dof) = gshape(spf, 2);
                        B(3, spf + dof) = gshape(spf, 2);
                        B(3, spf + 2 * dof) = gshape(spf, 1);
                        B(4, spf) = gshape(spf, 2);
                        B(4, spf + 2 * dof) = gshape(spf, 0);
                        B(5, spf) = gshape(spf, 1);
                        B(5, spf + dof) = gshape(spf, 0);
                    }
                }

            mfem::Vector temp(elstrain.Size()); temp = 0.0;
            B.Mult(eldofdisp,temp);
            double w = 1.0/(ir->GetNPoints());
            elstrain.Add(w,temp);
        }
    };

        void ElementStressStrain::ComputeElementStress(mfem::Vector &elstrain,  mfem::Coefficient &e, 
                                                        mfem::Coefficient &nu, mfem::Vector &elstress){

            int dim = el->GetDim();

            // Stiffness in Voigt notation. The stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            mfem::DenseMatrix C;  
            
            if (dim == 2)
            {
                C.SetSize(3, 3);
                elstress.SetSize(3);
            }
            else if (dim == 3)
            {
                C.SetSize(6, 6);
                elstress.SetSize(6);
            }
            
            C = 0.0;
            elstress = 0.0;

            elastic_mod = &e;
            poisson_ratio = &nu;

            const mfem::IntegrationRule *ir = GetIntegrationRule(*el, *Tr);
            if (ir == NULL)
            {
                int order = 2 * Tr->OrderGrad(el); // correct order?
                ir = &mfem::IntRules.Get(el->GetGeomType(), order);
            }

            mfem::real_t NU, E;
            mfem::DenseMatrix Ctemp;
            Ctemp.SetSize(C.NumRows(), C.NumCols());
            double w = 1.0/ir->GetNPoints();
        
            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                NU = E = 0;
                Ctemp = 0.0;
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                Tr->SetIntPoint(&ip);

                // The elastic constants are evaluated at each integration point.
                mfem::real_t NU_temp, E_temp;
                NU = poisson_ratio->Eval(*Tr, ip);
                E = elastic_mod->Eval(*Tr, ip); 

                if (dim == 2) 
                {
                    // Plane strain

                    // C(0, 0) = C(1, 1) = E / (1 - pow(NU,2));
                    // C(0, 1) = C(1, 0) = (E * NU) / (1 - pow(NU,2));
                    // C(2, 2) =  (E * (1 - NU)) / (2 * (1 - pow(NU,2)));

                    // Plane stress

                    C(0, 0) = C(1, 1) = (E /(1 - pow(NU,2)));
                    C(0, 1) = C(1, 0) = (E * NU /(1 - pow(NU,2)));
                    C(2, 2) =  (E * (1 - NU)/(2 * (1 - pow(NU,2))));

                    // 3D to 2D, but this is improper
                    
                    // C(0, 0) = C(1, 1) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    // C(0, 1) = C(1, 0) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    // C(2, 2) =  E / (2 * (1 + NU));
                }
                else if (dim == 3)
                {
                    Ctemp(0, 0) = Ctemp(1, 1) = Ctemp(2, 2) = (E * (1 - NU)) / ((1 - 2 * NU) * (1 + NU));
                    Ctemp(0, 1) = Ctemp(0, 2) = Ctemp(1, 0) = Ctemp(1, 2) = Ctemp(2, 0) = Ctemp(2, 1) = (E * NU) / ((1 - 2 * NU) * (1 + NU));
                    Ctemp(3, 3) = Ctemp(4, 4) = Ctemp(5, 5) =  E / (2 * (1 + NU));
                }

                C.Add(w, Ctemp);
            };
            
            C.Mult(elstrain, elstress);
        };

        void ElementStressStrain::ComputeElementStress(mfem::Vector &elstrain, mfem::MatrixCoefficient &Cmat, 
                                                    mfem::Vector &elstress){

           int dim = el->GetDim();

            // Stiffness in Voigt notation. The stiffness matrix has dimensions 3 x 3 in 2D and 6 x 6 in 3D.
            mfem::DenseMatrix C;  
            
            if (dim == 2)
                C.SetSize(3, 3);
            else if (dim == 3)
                C.SetSize(6, 6);

            elastic_constants = &Cmat;

            const mfem::IntegrationRule *ir = GetIntegrationRule(*el, *Tr);
            if (ir == NULL)
            {
                int order = 2 * Tr->OrderGrad(el); // correct order?
                ir = &mfem::IntRules.Get(el->GetGeomType(), order);
            }

            double w = 1.0/ir->GetNPoints();
            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                    const mfem::IntegrationPoint &ip = ir->IntPoint(i);

                    Tr->SetIntPoint(&ip);

                    mfem::DenseMatrix C_temp;
                    elastic_constants->Eval(C_temp, *Tr, ip);

                    C.Add(w, C_temp);
            }
            
            C.Mult(elstrain, elstress);
        };

        void ElementStressStrain::ComputeBoundaryElementArea(mfem::real_t &area){

            area = 0.0;
            mfem::real_t w = 0.0;
            int order = el->GetOrder();
            int geometry = el->GetGeomType();

            const mfem::IntegrationRule *ir = &mfem::IntRules.Get(geometry, order);

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const mfem::IntegrationPoint &ip = ir->IntPoint(i);
                Tr->SetIntPoint(&ip);
                w = ip.weight * Tr->Weight();                      // Quadrature weights
                area += w;
            }
        };

        void GlobalStressStrain::GlobalStrain(mfem::GridFunction &disp, mfem::GridFunction &strain){

            int numels = fespace->GetNE();
            int dim = mesh->Dimension();

            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;

            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;

            mfem::L2_FECollection* L2fec;
            mfem::FiniteElementSpace* L2fespace;

            L2fec = new mfem::L2_FECollection(0, dim);
            L2fespace = new mfem::FiniteElementSpace(mesh, L2fec, str_comp);

            strain = mfem::GridFunction(L2fespace);
            strain = 0.0;

            // Now start a loop that assembles strain vector element by element, and assembles the grid function.
            // The zeroth to numels index is \epsilon_{11}, then \epsilon_{22}, \epsilon_{33}, 
            // \epsilon_{23}, \epsilon_{13}, \epsilon_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = fespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = fespace->GetElementTransformation(elnum);

                mfem::Vector elstrain;
                ElementStressStrain Element(*fel, *Tr, elnum, fespace);
                Element.ComputeElementStrain(disp, elstrain);

                // elstrain.SetSize(6);
                // elstrain = 1;

                for (int comp = 0; comp < str_comp; comp++)
                    strain(elnum + (numels * comp)) = elstrain(comp);
    
            }
            
            delete L2fec;
            delete L2fespace;
        };

        void GlobalStressStrain::GlobalStrain(mfem::ParGridFunction &disp, mfem::ParGridFunction &strain){

            int numels = parfespace->GetNE();
            int dim = pmesh->Dimension();
            
            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;
            
            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;

            // Now start a loop that assembles strain vector element by element, and assembles the grid function.
            // The zeroth to numels index is \epsilon_{11}, then \epsilon_{22}, \epsilon_{33}, 
            // \epsilon_{23}, \epsilon_{13}, \epsilon_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = parfespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = parfespace->GetElementTransformation(elnum);

                mfem::Vector elstrain;
                ElementStressStrain Element(*fel, *Tr, elnum, parfespace);
                Element.ComputeElementStrain(disp, elstrain);

                // elstrain.SetSize(6);
                // elstrain = 1;

                for (int comp = 0; comp < str_comp; comp++)
                    strain(elnum + (numels * comp)) = elstrain(comp);
    
            }
        };


        void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain, mfem::Coefficient &e, 
                                            mfem::Coefficient &nu, mfem::GridFunction &stress){
            
            int numels = fespace->GetNE();
            int dim = mesh->Dimension();

            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;
            
            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;

            stress.SetSize(strain.Size());
            stress = 0.0;

            // Now start a loop that assembles stress vector element by element, and assembles the grid function.
            // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33}, 
            // \sigma_{23}, \sigma_{13}, \sigma_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = fespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = fespace->GetElementTransformation(elnum);
                mfem::Vector elstress(str_comp);
                mfem::Vector elstrain(str_comp);

                for (int comp = 0; comp < str_comp; comp++)
                    elstrain(comp) = strain(elnum + (numels * comp));
                
                ElementStressStrain Element(*fel, *Tr, elnum, fespace);
                Element.ComputeElementStress(elstrain, e, nu, elstress);

                for (int comp = 0; comp < str_comp; comp++)
                    stress(elnum + (numels * comp)) = elstress(comp);
    
            }
        };

        void GlobalStressStrain::GlobalStress(mfem::GridFunction &strain,
                                              mfem::MatrixCoefficient &Cmat, mfem::GridFunction &stress){
            
            int numels = fespace->GetNE();
            int dim = mesh->Dimension();

            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;
            
            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;

            stress.SetSize(strain.Size());
            stress = 0.0;

            // Now start a loop that assembles stress vector element by element, and assembles the grid function.
            // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33}, 
            // \sigma_{23}, \sigma_{13}, \sigma_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = fespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = fespace->GetElementTransformation(elnum);
                mfem::Vector elstress(str_comp);
                mfem::Vector elstrain(str_comp);

                for (int comp = 0; comp < str_comp; comp++)
                    elstrain(comp) = strain(elnum + (numels * comp));
                
                ElementStressStrain Element(*fel, *Tr, elnum, fespace);
                Element.ComputeElementStress(elstrain, Cmat, elstress);

                // elstrain.SetSize(6);
                // elstrain = 1;

                for (int comp = 0; comp < str_comp; comp++)
                    stress(elnum + (numels * comp)) = elstress(comp);
    
            }
        
        };

        void GlobalStressStrain::GlobalStress(mfem::ParGridFunction &strain, mfem::Coefficient &e, 
                                            mfem::Coefficient &nu, mfem::ParGridFunction &stress){
            
            int numels = parfespace->GetNE();
            int dim = pmesh->Dimension();
            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;
            
            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;
    
            // Now start a loop that assembles stress vector element by element, and assembles the grid function.
            // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33}, 
            // \sigma_{23}, \sigma_{13}, \sigma_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = parfespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = parfespace->GetElementTransformation(elnum);
                mfem::Vector elstress(str_comp);
                mfem::Vector elstrain(str_comp);

                for (int comp = 0; comp < str_comp; comp++)
                    elstrain(comp) = strain(elnum + (numels * comp));
                
                ElementStressStrain Element(*fel, *Tr, elnum, fespace);
                Element.ComputeElementStress(elstrain, e, nu, elstress);

                for (int comp = 0; comp < str_comp; comp++)
                    stress(elnum + (numels * comp)) = elstress(comp);
            }
        };

        void GlobalStressStrain::GlobalStress(mfem::ParGridFunction &strain,
                                              mfem::MatrixCoefficient &Cmat, mfem::ParGridFunction &stress){
            
            int numels = parfespace->GetNE();
            int dim = pmesh->Dimension();

            // Set up L2 FESpace to obtain one dof per element for each strain component. 
            // There are 3 strain components in 2D and 6 in 3D.
            int str_comp = 3;
            
            if (dim == 2)
                str_comp = 3;
            else if (dim == 3)
                str_comp = 6;

            stress.SetSize(strain.Size());
            stress = 0.0;

            // Now start a loop that assembles stress vector element by element, and assembles the grid function.
            // The zeroth to numels index is \sigma_{11}, then \sigma_{22}, \sigma_{33}, 
            // \sigma_{23}, \sigma_{13}, \sigma_{12}.

            for (int elnum = 0; elnum < numels; elnum++){

                const mfem::FiniteElement* fel = parfespace->GetFE(elnum);
                mfem::ElementTransformation* Tr = parfespace->GetElementTransformation(elnum);
                mfem::Vector elstress(str_comp);
                mfem::Vector elstrain(str_comp);

                for (int comp = 0; comp < str_comp; comp++)
                    elstrain(comp) = strain(elnum + (numels * comp));
                
                ElementStressStrain Element(*fel, *Tr, elnum, fespace);
                Element.ComputeElementStress(elstrain, Cmat, elstress);

                // elstrain.SetSize(6);
                // elstrain = 1;

                for (int comp = 0; comp < str_comp; comp++)
                    stress(elnum + (numels * comp)) = elstress(comp);
            }
        
        };

        void GlobalStressStrain::Global3DStrainComponents(mfem::GridFunction &strain, mfem::GridFunction &eps11, mfem::GridFunction &eps22,
                                    mfem::GridFunction &eps33, mfem::GridFunction &eps23, mfem::GridFunction &eps13, 
                                    mfem::GridFunction &eps12){

            mfem::L2_FECollection* L2fec;
            mfem::FiniteElementSpace* L2fespace;

            int dim = mesh->Dimension();

            L2fec = new mfem::L2_FECollection(0, dim);
            L2fespace = new mfem::FiniteElementSpace(mesh, L2fec);
            // auto L2fec = std::make_shared<mfem::L2_FECollection>(0, dim);
            // auto L2fespace = std::make_shared<mfem::FiniteElementSpace>(mesh, L2fec); // Should do something like this eventually.
            int NDofs = L2fespace->GetNDofs();

            mfem::GridFunction eps_temp(L2fespace);
            eps11 = eps22 = eps33 = eps23 = eps13 = eps12 = eps_temp;

            eps11 = eps22 = eps33 = eps23 = eps13 = eps12 = 0.0;

            for (int dof = 0; dof < NDofs; dof++){
                eps11(dof) = strain(dof);
                eps22(dof) = strain(dof + 1 * NDofs);
                eps33(dof) = strain(dof + 2 * NDofs);
                eps23(dof) = strain(dof + 3 * NDofs)/2; // Engineering shear strain to true shear strain.
                eps13(dof) = strain(dof + 4 * NDofs)/2; // Engineering shear strain to true shear strain.
                eps12(dof) = strain(dof + 5 * NDofs)/2; // Engineering shear strain to true shear strain.
            }

            delete L2fec;
            delete L2fespace;
            
        };

        void GlobalStressStrain::Global3DStressComponents(mfem::GridFunction &stress, mfem::GridFunction &sig11, mfem::GridFunction &sig22,
                                    mfem::GridFunction &sig33, mfem::GridFunction &sig23, mfem::GridFunction &sig13, 
                                    mfem::GridFunction &sig12){

        mfem::L2_FECollection* L2fec;
        mfem::FiniteElementSpace* L2fespace;

        int dim = mesh->Dimension();

        L2fec = new mfem::L2_FECollection(0, dim);
        L2fespace = new mfem::FiniteElementSpace(mesh, L2fec);
        int NDofs = L2fespace->GetNDofs();

        mfem::GridFunction sig_temp(L2fespace);
        sig11 = sig22 = sig33 = sig23 = sig13 = sig12 = sig_temp;

        sig11 = sig22 = sig33 = sig23 = sig13 = sig12 = 0.0;

        for (int dof = 0; dof < NDofs; dof++){
    
            sig11(dof) = stress(dof);
            sig22(dof) = stress(dof + 1 * NDofs);
            sig33(dof) = stress(dof + 2 * NDofs);
            sig23(dof) = stress(dof + 3 * NDofs);
            sig13(dof) = stress(dof + 4 * NDofs);
            sig12(dof) = stress(dof + 5 * NDofs);
        }
        }; 

        double GlobalStressStrain::ComputeBoundaryForce(mfem::GridFunction &stress, int &bdr_attribute, int &component)
        {
            double top_force = 0.0;
            
            int num_bdr_els = fespace->GetNBE();
            int num_els = fespace->GetNE();
            for (int i = 0; i < num_bdr_els; i++)
            {
                int bdr_attr = fespace->GetBdrAttribute(i);
                if (bdr_attr == bdr_attribute)
                {
                    mfem::FaceElementTransformations *ftr = fespace->GetMesh()->GetBdrFaceTransformations(i);
                    int volume_elem_index = ftr->Elem1No;
                    double element_sig = stress(volume_elem_index + (component - 1) * num_els);
                    double element_force = 0.0;
                    mfem::real_t element_area;
                    const mfem::FiniteElement &surf_el = *(fespace->GetBE(i));
                    mfem::ElementTransformation &el_trans = *(fespace->GetBdrElementTransformation(i));
                    mfemplus::ElementStressStrain ElementArea(surf_el, el_trans);
                    ElementArea.ComputeBoundaryElementArea(element_area);
                    element_force = element_area * element_sig;
                    top_force += element_force;
                }
            }
            return top_force;
   };

    };