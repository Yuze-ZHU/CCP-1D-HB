#include <iostream>

#include "ccp-1d.hpp"
#include "Const.hpp"
#include "../eigen-3.3.9/Eigen/Dense"
#include "../eigen-3.3.9/Eigen/Sparse"
void FVM::Solver::setupPoisson()
{
    using namespace Const;
    using Tools::toRowMajor;
    MatZeroEntries(poissonMat);
    // Solving A · x = b
    RowMajorMatrixXd A(numCells, numCells);

    A.setZero();
    // Construct matrix A
    for(label i = 0; i < numCells; ++i)
    {
        scalar invDxL   = 1.0 / faces[i].dist;      
        scalar invDxR   = 1.0 / faces[i + 1].dist;      
        A(i, i)         = eps0 * (invDxL + invDxR);          
        
        if (i > 0)
        {
            A(i, i - 1) =  -eps0 * invDxL;   
        }   

        if (i < numCells - 1)
        {
            A(i, i + 1) =  -eps0 * invDxR; 
        }                  
    }
             

    for (PetscInt iT = 0; iT < numT; ++iT)
    {
        PetscInt matBlockidx = iT;                               
        MatSetValuesBlocked(poissonMat,
                            1, &matBlockidx, 1, &matBlockidx,
                            A.data(), INSERT_VALUES);
    }

    MatAssemblyBegin(poissonMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(poissonMat,   MAT_FINAL_ASSEMBLY);
    
    // PetscViewer viewer;
    // PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);
    // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(poissonMat, viewer); 
    // PetscViewerPopFormat(viewer);  
}



void FVM::Solver::updatePhiFVM()
{
    using namespace Const;
    using Tools::toRowMajor;
    Eigen::VectorXd b(numCells);
    std::vector<PetscScalar> vecBlockVals(numCells);
    VecZeroEntries(phiRhsVec);
    for(label iT = 0; iT < numT; ++iT)
    {
        // Construct vector b
        b.setZero();
        for(label i = 0; i < numCells; ++i)
        {
            b(i) = e * (cells[i].Ni(iT) - cells[i].Ne(iT)) * cells[i].vol;

            if(i == 0)
            { 
                b(i) += eps0 * faces[i].PhiL(iT) / faces[i].dist;
            }

            if(i == numCells - 1)
            {
                b(i) += eps0 * faces[i + 1].PhiR(iT) / faces[i + 1].dist;
            }
        }
        
        std::memcpy(vecBlockVals.data(), b.data(), numCells * sizeof(PetscScalar));
        PetscInt vecBlockIdx = iT;
        VecSetValuesBlocked(phiRhsVec,1, &vecBlockIdx,          
                            vecBlockVals.data(), INSERT_VALUES);
    }
    VecAssemblyBegin(phiRhsVec);  
    VecAssemblyEnd  (phiRhsVec);
    // PetscViewer viewer;
    // PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);
    // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // VecView(phiRhsVec, viewer);   
    // PetscViewerPopFormat(viewer);

    // Create solver
    KSP phiSolver;
    // Create pre-conditioner
    PC phiPC;
    KSPCreate(PETSC_COMM_SELF, &phiSolver); 

    // Set coefficient matrix and pre-conditioner matrix                      
    KSPSetOperators(phiSolver, poissonMat, poissonMat); 
    
    // ---------- LU decomposition ------------
    KSPSetType(phiSolver, KSPPREONLY);
    KSPGetPC(phiSolver, &phiPC);
    PCSetType(phiPC, PCLU);

    // ---------------GMRES -------------------
    // KSPSetType(phiSolver, KSPGMRES);                                 
    // KSPGetPC(phiSolver, &phiPC);
    // PCSetType(phiPC, PCHYPRE);
    // PCHYPRESetType(phiPC,"boomeramg");   
    // KSPSetTolerances(phiSolver,
    //              1e-16,   // rtol
    //              1e-20,   // abstol
    //              1e5,     // dtol
    //              10000);   // max_iter



    KSPSetFromOptions(phiSolver);                                
    // Solving: x = A⁻¹·b
    KSPSolve(phiSolver, phiRhsVec, phiVec); 


    // Assign phi value to cell
    PetscScalar *phiArr = nullptr;
    VecGetArray(phiVec, &phiArr);                    
    for (PetscInt iT = 0; iT < numT; ++iT)
    {
        PetscInt offset = iT * numCells;       

        for (PetscInt i = 0; i < numCells; ++i)
        {
            cells[i].Phi(iT) = phiArr[offset + i];    
        }

        for (label i = 0; i < numCells; ++i)
        {
            if (i == 0)
            {
                cells[i].Ec(iT) = -(cells[i].Phi(iT) - faces[i].PhiL(iT)) / faces[i].dist;
            }
            else if (i == numCells - 1)
            {
                cells[i].Ec(iT) = -(faces[i + 1].PhiR(iT) - cells[i].Phi(iT)) / faces[i + 1].dist;
            }
            else
            {
                cells[i].Ec(iT) = -(cells[i + 1].Phi(iT) - cells[i - 1].Phi(iT))
                                / (faces[i].dist + faces[i +1].dist);
            }
        }
        // Ghost cells at the boundaries
        cells[numCells].Ec(iT)     = cells[numCells - 1].Ec(iT);
        cells[numCells + 1].Ec(iT) = cells[0].Ec(iT);
    }
    VecRestoreArray(phiVec, &phiArr);         
    KSPDestroy(&phiSolver);   
}

