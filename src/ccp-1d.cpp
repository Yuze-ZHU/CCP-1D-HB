// -------------------------------------------------------//
//  _______  _______  _______          __    ______       //
// (  ____ \(  ____ \(  ____ )        /  \  (  __  \      //
// | (    \/| (    \/| (    )|        \/) ) | (  \  )     //
// | |      | |      | (____)| _____    | | | |   ) |     //
// | |      | |      |  _____)(_____)   | | | |   | |     //
// | |      | |      | (                | | | |   ) |     //
// | (____/\| (____/\| )              __) (_| (__/  )     //
// (_______/(_______/|/               \____/(______/      //
//                                                        //
// -------------------------------------------------------//
//                                                        //
// Description:                                           //
// Contains Const, Solver, Tools                          //
//                                                        //
// -------------------------------------------------------//

#include <iostream>
#include <cmath>
#include <unordered_map>
#include <chrono>

// Output part
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include "ccp-1d.hpp"
#include "Const.hpp"
#include "def.hpp"
#include <cmath>

// ------------------------------------------------------ //
// Description:                                           //
// Tools for basic kinetic theory                         //
// ------------------------------------------------------ //

using std::vector;

namespace Tools
{

    scalar abs(const scalar num)
    {
        return num > 0 ? num : -num;
    }

    tensor1 sign(const tensor1 &array)
    {
        tensor1 sign = array;
        sign[sign >= 0.0] = 1.0;
        sign[sign < 0.0] = -1.0;

        return sign;
    }

    scalar sign(const scalar &value)
    {
        return value > 0.0 ? 1.0 : -1.0;
    }

    scalar max(const scalar a, const scalar b)
    {
        return a > b ? a : b;
    }

    scalar min(const scalar a, const scalar b)
    {
        return a < b ? a : b;
    }

    scalar KToeV(const scalar K)
    {
        return K * 8.6174e-5;
    }

    scalar eVtoK(const scalar eV)
    {
        return eV * 1.1604e4;
    }

    scalar TeToEe(const scalar Ne, const scalar Te)
    {
        return 1.5 * Ne * Const::kB * Te;
    }


    scalar TeToEeND(const scalar Ne, const scalar Te)
    {
        return Ne * Te;
    }

    scalar EeToTeND(const scalar Ne, const scalar Ee)
    {
        return Ee / Ne;
    }

    scalar EeToTe(const scalar Ne, const scalar Ee)
    {
        return 0.6666667 * Ee / ( Ne*Const::kB);
    }

    std::string toString(double value, int precision) 
    {
        std::ostringstream out;
        out << std::fixed 
            << std::setprecision(precision) << value;
        return out.str();
    }

    std::string toStringFill0(double value)
    {
        //- The Func setw() controls the digits before the decimal point
        //- The Func fixed() with setprecision() 
        //- controls the digits before the decimal point
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(4) 
           << std::fixed << std::setprecision(2) << value;

        std::string result = ss.str();
        return result;
    }

    std::vector<PetscScalar> toRowMajor(const Eigen::MatrixXd& M)
    {
        std::vector<PetscScalar> buf(M.rows() * M.cols());

        for (PetscInt r = 0; r < M.rows(); ++r)
            for (PetscInt c = 0; c < M.cols(); ++c)
                buf[r * M.cols() + c] = static_cast<PetscScalar>( M(r, c) ); 

        return buf; 
    }
}

// ------------------------------------------------------- //
// Description:                                            //
// Constructor and deconstructor for cell class            //
// ------------------------------------------------------- //

FVM::ChemSet::ChemSet(const std::string csvName)
{
    using std::ifstream;
    using std::string;
    using std::stringstream;
    using std::vector;

    ifstream chemCSV(csvName);

    if ( chemCSV.is_open() )
    {
        std::string word;
        if (chemCSV >> word) {
            std::cout << "Successfully Load Chemical Data: " 
            << word << std::endl;
        } else {
            std::cout << "Cannot Load Chemical Data." << std::endl;
        }
    }
    else
    {
        std::cout << "Cannot open the .csv file." << std::endl; 
    }

    string datas;
    vector<string> headers;

    if (getline(chemCSV, datas))
    {
        stringstream ss(datas);
        string word;
        while (getline(ss, word, ','))
        {
            headers.push_back(word);
        }
    }

    while (getline(chemCSV, datas))
    {
        stringstream ss(datas);
        string data;
        vector<string> row;

        while (getline(ss, data, ','))
        {
            row.push_back(data);
        }

        if (row.size() == 2)
        {
            ChemRate dp;
            dp.energy = stod(row[0]);
            dp.cRate = stod(row[1]);
            chemSets.push_back(dp);
        }
    }

    chemCSV.close();
}

FVM::ChemSet::~ChemSet()
{}

scalar FVM::ChemSet::interpolateChem(const scalar E)
{
    size_t low = 0;
    size_t high = chemSets.size() - 1;

    // Bindary search to find the interval where the energy E lies
    while (low < high)
    {
        size_t mid = low + (high - low) * 0.5;

        if (chemSets[mid].energy < E)
        {
            low = mid + 1;
        }
        else
        {
            high = mid;
        }
    }

    if (low == 0)
    {
        return chemSets[0].cRate;
    }
    else if 
    (
        low == chemSets.size() || 
        chemSets[low].energy == E
    )
    {
        return chemSets[low].cRate;
    }
    else
    {
        const double E0 = chemSets[low - 1].energy;
        const double E1 = chemSets[low].energy;
        const double k0 = chemSets[low - 1].cRate;
        const double k1 = chemSets[low].cRate;

        // Linear interpolation
        return k0 + (k1 - k0) * (E - E0) / (E1 - E0);
    }

}


FVM::Cell::Cell(const label i)
    : vol(0.0), cellID(i), dt(0.0), dtau(0.0), 
    sumNe(0.0), sumNi(0.0), sumEe(0.0), sumTe(0.0), sumPhi(0.0), sumE(0.0),
    NeA(0.0), NiA(0.0), EeA(0.0), TeA(0.0), PhiA(0.0), EA(0.0),

    // Vectors initialization
    // Rate coefficient
    kl(Eigen::VectorXd::Zero(Const::numT)),
    // Conservative variables
    Ne(Const::numT), Ni(Const::numT), Ee(Const::numT), Phi(Const::numT),

    // Primitive variables
    Te(Const::numT), Ec(Eigen::VectorXd::Zero(Const::numT)),

    // Other variables
    Je(Eigen::VectorXd::Zero(Const::numT)), 
    Ji(Eigen::VectorXd::Zero(Const::numT)), 
    cS(Eigen::VectorXd::Zero(Const::numT)),

    // Conservative variables at the previous time step
    NeOld(Const::numT), NiOld(Const::numT), EeOld(Const::numT), PhiOld(Const::numT), 

    //- Primitive variables at the previous time step
    TeOld(Const::numT),    
    
    //- Residual flux of conservative variables 
    //- (including flux and source term, updating the conservative variables)
    ResFluxNe(Eigen::VectorXd::Zero(Const::numT)),
    ResFluxNi(Eigen::VectorXd::Zero(Const::numT)),
    ResFluxEe(Eigen::VectorXd::Zero(Const::numT)),
    ResFluxPhi(Eigen::VectorXd::Zero(Const::numT)),
    
    //- Slope of conservative variables
    sNe(Eigen::VectorXd::Zero(Const::numT)),
    sNi(Eigen::VectorXd::Zero(Const::numT)),
    sEe(Eigen::VectorXd::Zero(Const::numT)),
    sPhi(Eigen::VectorXd::Zero(Const::numT)),

    //- Stiffness matrix
    massDiagnal(RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS))

{
    // Matrix for the chemical reaction
    jacobianCs.resize(Const::numT);
    for (label iT = 0; iT < Const::numT; ++iT)
    {
        Ne(iT)     = Const::Ne0; 
        Ni(iT)     = Const::Ni0;
        Te(iT)     = Const::Te0;
        Ee(iT)     = Tools::TeToEeND(Ne(iT), Te(iT));   
        Phi(iT)    = Const::Phi0; 

        NeOld(iT)  = Ne(iT);
        NiOld(iT)  = Ni(iT);
        TeOld(iT)  = Te(iT);
        EeOld(iT)  = Ee(iT);
        PhiOld(iT) = Phi(iT);

        jacobianCs[iT] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
    }
}


FVM::Cell::~Cell()
{

}



// ------------------------------------------------------- //
// Description:                                            //
// Constructor and deconstructor for face class            //
// ------------------------------------------------------- //

FVM::Face::Face(const Cell &Lcell, const Cell &Rcell)
    : cellL(Lcell), cellR(Rcell), faceID(Rcell.cellID),
    cellIDL(Lcell.cellID), cellIDR(Rcell.cellID), dist(0.0),
    distL(0.0), distR(0.0),

    //- Vectors initialization
    //- Conservative and primitive variables at left and right sides
    NeL(Eigen::VectorXd::Zero(Const::numT)),
    NeR(Eigen::VectorXd::Zero(Const::numT)),
    NiL(Eigen::VectorXd::Zero(Const::numT)),
    NiR(Eigen::VectorXd::Zero(Const::numT)),
    EeL(Eigen::VectorXd::Zero(Const::numT)),
    EeR(Eigen::VectorXd::Zero(Const::numT)),
    TeL(Eigen::VectorXd::Zero(Const::numT)),
    TeR(Eigen::VectorXd::Zero(Const::numT)),
    PhiL(Eigen::VectorXd::Zero(Const::numT)),
    PhiR(Eigen::VectorXd::Zero(Const::numT)),
    Ef(Eigen::VectorXd::Zero(Const::numT)),

    //- Flux for number density and electron energy
    fluxNe(Eigen::VectorXd::Zero(Const::numT)),
    fluxNi(Eigen::VectorXd::Zero(Const::numT)),
    fluxEe(Eigen::VectorXd::Zero(Const::numT)),
    fluxEeJoule(Eigen::VectorXd::Zero(Const::numT)),
    fluxPhi(Eigen::VectorXd::Zero(Const::numT))
{
    
    if (faceID == 0)
    {
        cellIDL = -1;
    }
    else if (faceID == Const::numCells)
    {
        cellIDR = -1;
    }

    //- Flux Jacobian
    jacobianFluxL.resize(Const::numT);
    jacobianFluxR.resize(Const::numT);
    jacobianFluxLNeg.resize(Const::numT);
    jacobianFluxRNeg.resize(Const::numT);
    for (label iT = 0; iT < Const::numT; ++iT)
    {

        jacobianFluxL[iT]    = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
        jacobianFluxR[iT]    = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
        jacobianFluxLNeg[iT] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
        jacobianFluxRNeg[iT] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
    }
}




FVM::Face::~Face()
{

}

scalar FVM::Face::getDiff(const scalar pLcell, const scalar pRcell)
{
    // Diffusion gets the gradient at the cell center
    return (pRcell - pLcell) / dist;
}

scalar FVM::Face::getUpwind(const scalar pLcell, const scalar pRcell, const label iT)
{
    // Upwind uses variables at the cell interface
    return cellL.Phi(iT) > cellR.Phi(iT) ? pRcell : pLcell;
}

scalar FVM::Face::getDownwind(const scalar pLcell, const scalar pRcell, const label iT)
{
    // Upwind uses variables at the cell interface
    return cellL.Phi(iT) > cellR.Phi(iT) ? pLcell : pRcell;
}

void FVM::Face::interpolate()
{
    using namespace Const;
    using namespace Tools;
    
    for (label iT = 0; iT < numT; ++iT)
    {
        if(faceID == 0)
        {
            NeR(iT)  = cellR.Ne(iT)  - 0.5 * distR * cellR.sNe(iT);
            NiR(iT)  = cellR.Ni(iT)  - 0.5 * distR * cellR.sNi(iT);
            EeR(iT)  = cellR.Ee(iT)  - 0.5 * distR * cellR.sEe(iT);
            PhiR(iT) = cellR.Phi(iT) - 0.5 * distR * cellR.sPhi(iT);
            TeR(iT)  = EeToTeND(NeR(iT), EeR(iT));

            // For the boundary cells, distL and slopes are zero
            NeL(iT)  = cellL.Ne(iT)  + 0.5 * distL * cellL.sNe(iT);
            NiL(iT)  = cellL.Ni(iT)  + 0.5 * distL * cellL.sNi(iT);
            EeL(iT)  = cellL.Ee(iT)  + 0.5 * distL * cellL.sEe(iT);
            PhiL(iT) = cellL.Phi(iT) + 0.5 * distL * cellL.sPhi(iT);
            // Zero gradient for electron energy at the boundary
            TeL(iT)  = TeR(iT);
        }
        else if(faceID == numCells)
        {
            NeL(iT)  = cellL.Ne(iT)  + 0.5 * distL * cellL.sNe(iT);
            NiL(iT)  = cellL.Ni(iT)  + 0.5 * distL * cellL.sNi(iT);
            EeL(iT)  = cellL.Ee(iT)  + 0.5 * distL * cellL.sEe(iT);
            PhiL(iT) = cellL.Phi(iT) + 0.5 * distL * cellL.sPhi(iT);
            TeL(iT)  = EeToTeND(NeL(iT), EeL(iT));
            // For the boundary cells, distR and slopes are zero
            NeR(iT)  = cellR.Ne(iT)  - 0.5 * distR * cellR.sNe(iT);
            NiR(iT)  = cellR.Ni(iT)  - 0.5 * distR * cellR.sNi(iT);
            EeR(iT)  = cellR.Ee(iT)  - 0.5 * distR * cellR.sEe(iT);
            PhiR(iT) = cellR.Phi(iT) - 0.5 * distR * cellR.sPhi(iT);
            // Zero gradient for electron energy at the boundary
            TeR(iT)  = TeL(iT);
        }
        else
        {
            NeL(iT)  = cellL.Ne(iT)  + 0.5 * distL * cellL.sNe(iT);
            NiL(iT)  = cellL.Ni(iT)  + 0.5 * distL * cellL.sNi(iT);
            EeL(iT)  = cellL.Ee(iT)  + 0.5 * distL * cellL.sEe(iT);
            PhiL(iT) = cellL.Phi(iT) + 0.5 * distL * cellL.sPhi(iT);
            TeL(iT)  = EeToTeND(NeL(iT), EeL(iT));

            NeR(iT)  = cellR.Ne(iT)  - 0.5 * distR * cellR.sNe(iT);
            NiR(iT)  = cellR.Ni(iT)  - 0.5 * distR * cellR.sNi(iT);
            EeR(iT)  = cellR.Ee(iT)  - 0.5 * distR * cellR.sEe(iT);
            PhiR(iT) = cellR.Phi(iT) - 0.5 * distR * cellR.sPhi(iT);
            TeR(iT)  = EeToTeND(NeR(iT), EeR(iT));
        }
    }
}

FVM::Solver::Solver()
    : 
    step(0), runTime(0.0), runCyc(0.0), dtminGlobal(0.0), 
    isWarning(false),outputDir("."), chemSets(Const::chemCSV),

    // Residual vectors for output
    resNe(Eigen::VectorXd::Zero(Const::numT)),
    resNi(Eigen::VectorXd::Zero(Const::numT)),
    resEe(Eigen::VectorXd::Zero(Const::numT)),
    resTe(Eigen::VectorXd::Zero(Const::numT)),

    //- Time spectral source matrix, E
    harmMat(RowMajorMatrixXd::Zero(Const::numT, Const::numT)),

    //- DFT matrix, D
    dftMat(RowMajorMatrixXd::Zero(Const::numT, Const::numT)),

    //- Inverse of DFT matrix, D^-1
    dftMatInv(RowMajorMatrixXd::Zero(Const::numT, Const::numT)),

    //- Time instants
    harmTime(Eigen::VectorXd::Zero(Const::numT)),

    //- Angular frequencies
    angFreq(Eigen::VectorXd::Zero(Tools::max(Const::numH, 1))),

    //- Harmonic balance source Jacobian
    jacobianHBs(RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS))    

{
    std::cout << " ------ Starting constructing solver ------ " << std::endl;
    std::cout << " ------ Constucting cells ------" << std::endl;
    for (label i = 0; i < Const::numCells + 2; i++)
    {
        cells.emplace_back(i);
    }
    std::cout << " ------ Constructing cells finished ------" << std::endl;
    std::cout << " ------ Constructing faces --------" << std::endl;
    for (label i = 0; i < Const::numNodes; i++)
    {
        if (i == 0) // cell at left boundary
        {
            faces.emplace_back(cells.at(Const::numCells + 1), cells.at(i));
        }
        else // cells at the interior faces and right boundary
        {
            faces.emplace_back(cells.at(i - 1), cells.at(i));
        }
    }
    std::cout << " ------ Constructing faces finished ------" << std::endl;

    if(Const::implicitScheme == ImplicitScheme::NO ||
       Const::implicitScheme == ImplicitScheme::PCI1)
    {
        // Total number of global degrees of freedom(DoF)
        PetscInt numPhiDoF = Const::numCells * Const::numT;
        
        // (1):  Create matrix coefficient matrix A for the discretized Poisson's equation
        // based on PETSc in sequential (non-parallel) mode
        MatCreate(PETSC_COMM_SELF, &poissonMat);

        // Set matrix poissonMat global size: numPhiDoF × numPhiDoF 
        // Only one process owns all rows/cols
        MatSetSizes(poissonMat, PETSC_DECIDE, PETSC_DECIDE, numPhiDoF, numPhiDoF);

        // Set block size: each nonzero block is a numCells * numCells submatrix
        MatSetBlockSize(poissonMat, Const::numCells);

        // Use sequential block AIJ (BAIJ) format (serial block-sparse matrix)
        MatSetType(poissonMat, MATSEQBAIJ);

        // Preallocate memory: assume maximum 1 block entries per row (
        MatSeqBAIJSetPreallocation(poissonMat, Const::numCells, 1, nullptr);

        // Setup matrix
        MatSetUp(poissonMat);   

        // (2):  Create Solution vector x
        VecCreate(PETSC_COMM_SELF, &phiVec);
        VecSetSizes(phiVec, PETSC_DECIDE, numPhiDoF); 
        VecSetBlockSize(phiVec, Const::numCells);   
        VecSetFromOptions(phiVec);
        VecSetUp(phiVec);

        // (3):  Create RHS vector b
        VecDuplicate(phiVec, &phiRhsVec);
        std::cout << " ------ Allocating memeory for solving Poisson's equation Finished ------" << std::endl;
    }

    
    if ((Const::implicitScheme == ImplicitScheme::PCI1)
       ||(Const::implicitScheme == ImplicitScheme::FCI))
       {
            std::cout << " ------ Constructing matrices and vectors for implicit scheme --------" << std::endl;
            // Matrix creation
            MatCreate(PETSC_COMM_SELF, &jacobianAllGlobal);

            MatSetSizes(jacobianAllGlobal, PETSC_DECIDE, PETSC_DECIDE, Const::blockSizeL, Const::blockSizeL);

            MatSetBlockSize(jacobianAllGlobal, Const::blockSizeS);

            MatSetType(jacobianAllGlobal, MATSEQBAIJ);

            MatSeqBAIJSetPreallocation(jacobianAllGlobal, Const::blockSizeS, Const::numT + 2, nullptr);

            MatSetUp(jacobianAllGlobal); 


            // Vector creation
            VecCreate(PETSC_COMM_SELF, &incrementGlobal);

            VecSetSizes(incrementGlobal, PETSC_DECIDE, Const::blockSizeL);

            VecSetBlockSize(incrementGlobal, Const::blockSizeS);

            VecSetFromOptions(incrementGlobal);
        
            VecSetUp(incrementGlobal);

            VecDuplicate(incrementGlobal, &rhsGlobal);

            std::cout << " ------ Constructing matrices and vectors for implicit scheme finished --------"
                      << std::endl;
       }
    
    jacobianFluxJouleLCell.resize(Const::numT);
    jacobianFluxJouleCCell.resize(Const::numT);
    jacobianFluxJouleRCell.resize(Const::numT);
    for (label iT = 0; iT < Const::numT; ++iT) 
    {
        jacobianFluxJouleLCell[iT].resize(Const::numCells);
        jacobianFluxJouleCCell[iT].resize(Const::numCells);
        jacobianFluxJouleRCell[iT].resize(Const::numCells);       
        for (label i = 0; i < Const::numCells; ++i) 
        {
            jacobianFluxJouleLCell[iT][i] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
            jacobianFluxJouleCCell[iT][i] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);
            jacobianFluxJouleRCell[iT][i] = RowMajorMatrixXd::Zero(Const::blockSizeS, Const::blockSizeS);        
        }
    }

    std::cout << " ------ Finish Construct Solver" << std::endl;
}

FVM::Solver::~Solver()
{

}

inline void FVM::Solver::addValuesBlock(Mat &mat, PetscInt row, PetscInt col,
                         const Eigen::Matrix<scalar, 4, 4, Eigen::RowMajor> &block)
{
    MatSetValuesBlocked(mat, 1, &row, 1, &col, block.data(), ADD_VALUES);
}


void FVM::Solver::initlize()
{
    using namespace Const;

    if (initOption != InitOption::SCRATCH) {
        // std::cout << "Read quasi-steady solution from file." << std::endl;
        // readQuasiSteadySolution(initFile);
    }

}

void FVM::Solver::initHarmonicMat()
{
    using namespace Const;
    const double freq = fRF;
    const double T    = period;
    double dtheta;

    RowMajorMatrixXd D     = RowMajorMatrixXd::Zero(numT, numT);
    RowMajorMatrixXd Dt    = RowMajorMatrixXd::Zero(numT, numT);
    RowMajorMatrixXd D_inv = RowMajorMatrixXd::Zero(numT, numT);


    // Initialize the harmonic time vector, angular frequency vector
    harmTime.setZero();
    angFreq.setZero();

    // Initialize the harmonic matrix, DFT matrix, and IDFT matrix
    harmMat.setZero();
    dftMat.setZero();
    dftMatInv.setZero();

    std::cout << "Time instants used in the harmonic balance analysis are:" << std::endl;
    for (label nT = 0; nT < numT; ++nT)
    {
        harmTime(nT) = nT * T / numT;
        std::cout << "t(" << nT << ") = " << harmTime(nT) << std::endl;
    }

    if (analysisMode == AnalysisMode::HB)
    {
        std::cout << "Angular frequencies used in the harmonic balance analysis are:\n";
        for (label nh = 0; nh < numH; ++nh)
        {
            angFreq(nh) = 2.0 * pi * (nh + 1) * freq;
            std::cout << "  angFreq(" << nh << ") = " << angFreq(nh) << '\n';
        }
    }
    else 
    {
        angFreq(0) = 0.0;
        std::cout << "Harmonic balance disabled. angFreq(0) = 0.0\n";
    }


    for (label nT = 0; nT < numT; ++nT)
    {
        D(nT, 0)  = 1.0;
        Dt(nT, 0) = 0.0;
        for (label nh = 0; nh < numH; ++nh)
        {   
            dtheta             = angFreq(nh) * harmTime(nT);
            D(nT, 2 * nh + 1)  = sin(dtheta);
            D(nT, 2 *nh + 2)   = cos(dtheta);
            // Nondimensionlization for angFreq
            Dt(nT, 2 * nh + 1) = angFreq(nh) * cos(dtheta);
            Dt(nT, 2 * nh + 2) = -angFreq(nh) * sin(dtheta);
        }
    }

    std::cout << "\nInverse Discrete Fourier transform operator D:" << std::endl;
    for (label i = 0; i < numT; ++i) {
        for (label j = 0; j < numT; ++j) {
            std::cout << std::setw(20) 
                      << std::fixed 
                      << std::setprecision(8) 
                      << D(i, j) << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "\nFirst-order derivative opeartot Dt:" << std::endl;
    for (label i = 0; i < numT; ++i) {
        for (label j = 0; j < numT; ++j) {
            std::cout << std::setw(20)
                      << std::fixed 
                      << std::setprecision(8) 
                      << Dt(i, j) << '\t';
        }
        std::cout << std::endl;
    }

    D_inv = D.inverse();

    std::cout << "\nDiscrete Fourier transform operator D^-1:" << std::endl;
    for (label i = 0; i < numT; ++i) 
    {
        for (label j = 0; j < numT; ++j) 
        {
            std::cout << std::setw(20) 
                      << std::fixed 
                      << std::setprecision(8) 
                      << D_inv(i, j) << '\t';
        }
        std::cout << std::endl;
    }

    harmMat = Dt * D_inv;
    
    std::cout << "\nTime spectral source operator E:" << std::endl;
    for (label i = 0; i < numT; ++i) {
        for (label j = 0; j < numT; ++j) {
            std::cout << std::setw(20) 
                      << std::fixed 
                      << std::setprecision(8) 
                      << harmMat(i, j) << '\t';
        }
        std::cout << std::endl;
    }

    dftMat    = D_inv; 
    dftMatInv = D;
}

void FVM::Solver::initRK()
{
    using namespace Const;
    for (label i = 0; i < numCells; ++i)
    {
        cells[i].NeOld   = cells[i].Ne;
        cells[i].NiOld   = cells[i].Ni;
        cells[i].EeOld   = cells[i].Ee;
        cells[i].PhiOld  = cells[i].Phi;
        cells[i].TeOld   = cells[i].Te; 
    }
}

void FVM::Solver::getDt(bool& isWritingStepCal, bool& isResWritingStepCal, bool& isPrintStepCal)
{
    using namespace Const;
    using Tools::min;

    dtminGlobal = 1e150;

    for (label i = 0; i < numCells; ++i)
    {
        scalar dtminLocal   = 1e150;
        const scalar dtDiff = 0.25 * cells[i].vol * cells[i].vol / De;
        for(label nT = 0; nT < numT; ++nT)
        {
            const scalar dtConv = cells[i].vol / (muE * abs(cells[i].Ec[nT]) + 1e-120 );
            scalar dtmin1       = min(dtDiff, dtConv);
            dtminLocal          = min(dtminLocal, dtmin1);   
        }
        // Modified local time step based on the biggest harmonic frequency
        dtminLocal  = CFL * (1. / (1. / dtminLocal + 2.0 * pi * fRF * numH));
        cells[i].dt = dtminLocal; 
        dtminGlobal = min(cells[i].dt, dtminGlobal);
    }
    
    if (timeStepType == TimeStepType::GLOBAL)
    {
        for (label i = 0; i < numCells; ++i)
        {
            cells[i].dt = dtminGlobal;
        }
    }
    
    // Convert cycle to step approximately
    // Calculated only once
    if (writeMode == WriteMode::CYC){
        if(!isWritingStepCal)
        {
            writeStep = int(period / dtminGlobal * writeCyc);
            std::cout << "writeStep is set to " << writeStep << std::endl;
            isWritingStepCal = true;
        }
    }

    if (resWriteMode == ResWriteMode::CYC){
        if(!isResWritingStepCal)
        {
            resWriteStep = int(period / dtminGlobal * resWriteCyc);
            std::cout << "resWriteStep sis set to " << resWriteStep << std::endl;
            isResWritingStepCal = true;
        }
    }

    if (printMode == PrintMode::CYC){
        if(!isPrintStepCal)
        {
            printStep = int(period / dtminGlobal * printCyc);
            std::cout << "printStep sis set to " << printStep << std::endl;
            isPrintStepCal = true;
        }
    }
}

void FVM::Solver::getDtau()
{
    // Physical time and Pseudo time share the same notation dt;
    using namespace Const;
    using Tools::min;

    dtminGlobal = 1e150;

    for (label i = 0; i < numCells; ++i)
    {
        scalar dtminLocal   = 1e150;
        const scalar dtDiff = 0.25 * cells[i].vol * cells[i].vol / De;
        for(label nT = 0; nT < numT; ++nT)
        {
            const scalar dtConv = cells[i].vol / (muE * abs(cells[i].Ec[nT]) + 1e-120 );
            scalar dtmin1       = min(dtDiff, dtConv);
            dtminLocal          = min(dtminLocal, dtmin1);   
        }
        // Modified local time step based on the biggest harmonic frequency
        dtminLocal  = CFL * (1. / (1. / dtminLocal + 2.0 * pi * fRF * numH));
        cells[i].dt = dtminLocal; 
        dtminGlobal = min(cells[i].dt, dtminGlobal);
    }
    
    if (timeStepType == TimeStepType::GLOBAL)
    {
        for (label i = 0; i < numCells; ++i)
        {
            cells[i].dt = dtminGlobal;
        }
    }
}

void FVM::Solver::getMassDiagnal()
{
    using namespace Const;
    for (label i = 0; i < numCells; ++i)
    {
        for (label bs = 0; bs < blockSizeS; ++bs)
        {
            cells[i].massDiagnal(bs, bs) = cells[i].vol / cells[i].dt;
        }
    }
}

void FVM::Solver::assembleGlobalVecRHS(const label iRK)
{
    using namespace Const;
    using namespace Tools;
    VecSet(rhsGlobal, 0.0);
    scalar TeV;

    for (label iT = 0; iT < numT; ++iT)
    { 
        for (label i = 0; i < numCells; ++i)
        {
            PetscScalar rhsiTCell[4];
            PetscInt blockID = iT * numCells + i;

            TeV = KToeV(EeToTe(cells[i].Ne(iT) * nRef, cells[i].Ee(iT) * EeRef));
            cells[i].kl(iT) = chemSets.interpolateChem(TeV) / klRef;

            cells[i].cS(iT)          = cells[i].kl(iT) * cells[i].Ne(iT) * N;
            cells[i].Je(iT)          = (faces[i].fluxNe(iT) + faces[i + 1].fluxNe(iT)) / 2.0;
            cells[i].Ji(iT)          = (faces[i].fluxNi(iT) + faces[i + 1].fluxNi(iT)) / 2.0;

            //- Flux terms
            cells[i].ResFluxNe(iT)   = alpha[iRK]*(faces[i].fluxNe(iT) - faces[i + 1].fluxNe(iT));
            cells[i].ResFluxNi(iT)   = alpha[iRK]*(faces[i].fluxNi(iT) - faces[i + 1].fluxNi(iT));
            cells[i].ResFluxEe(iT)   = alpha[iRK]*(faces[i].fluxEe(iT) - faces[i + 1].fluxEe(iT));
            cells[i].ResFluxEe(iT)  += alpha[iRK]*(faces[i].fluxEeJoule(iT) - faces[i + 1].fluxEeJoule(iT)) 
                                     * cells[i].Phi(iT);
            cells[i].ResFluxPhi(iT)  = alpha[iRK] * (faces[i].fluxPhi(iT) - faces[i + 1].fluxPhi(iT));

            //- Source terms
            cells[i].ResFluxNe(iT)  += alpha[iRK] * cells[i].kl(iT) * N * cells[i].Ne(iT) * cells[i].vol; 
            cells[i].ResFluxNi(iT)  += alpha[iRK] * cells[i].kl(iT) * N * cells[i].Ne(iT) * cells[i].vol;
            cells[i].ResFluxEe(iT)  -= alpha[iRK] * cells[i].kl(iT) * N * cells[i].Ne(iT) * Hl * cells[i].vol;
            cells[i].ResFluxPhi(iT) += alpha[iRK] * (cells[i].Ne(iT) - cells[i].Ni(iT))* cells[i].vol * e;

            rhsiTCell[0] = cells[i].ResFluxNe(iT);
            rhsiTCell[1] = cells[i].ResFluxNi(iT);
            rhsiTCell[2] = cells[i].ResFluxEe(iT);
            rhsiTCell[3] = cells[i].ResFluxPhi(iT);

            VecSetValuesBlocked(rhsGlobal, 1, &blockID, rhsiTCell, INSERT_VALUES);

        }
    }
    VecAssemblyBegin(rhsGlobal);
    VecAssemblyEnd(rhsGlobal);
}

void FVM::Solver::assembleLocalFluxJacobian()
{
    using namespace Const;
    for (auto &face : faces)
    {
        if (face.faceID == 0)
        {
            face.getFluxJacobianLeftBCFCI();
            face.getFluxJacobianNeg();
        }
        else if (face.faceID == numCells)
        {
            face.getFluxJacobianRightBCFCI();
            face.getFluxJacobianNeg();
        }
        else
        {
            face.getFluxJacobianFCI();
            face.getFluxJacobianNeg();
        }
    }    
}

void FVM::Solver::assembleGlobalJacobian()
{
    using namespace Const;
    using namespace Tools;
    DiagnalMatrixXd4 identity4 = DiagnalMatrixXd4::Identity();
    MatZeroEntries(jacobianAllGlobal);
    for (label iT = 0; iT < numT; ++iT)
    {
        // --- 1.  Flux Jacobian ---
        // --- 1.1 Flux Jacobian from the interior faces ---
        for (label f = 1; f < numCells; ++f)
        {
            const auto &face = faces[f];

            PetscInt idL = iT * numCells + face.cellIDL;
            PetscInt idR = iT * numCells + face.cellIDR;

            addValuesBlock(jacobianAllGlobal, idL, idL, face.jacobianFluxL[iT]);
            addValuesBlock(jacobianAllGlobal, idL, idR, face.jacobianFluxR[iT]);
            addValuesBlock(jacobianAllGlobal, idR, idL, face.jacobianFluxLNeg[iT]);
            addValuesBlock(jacobianAllGlobal, idR, idR, face.jacobianFluxRNeg[iT]);
        }
        // --- 1.2 Flux Jacobian from the left boundary faces ---       
        const auto &faceLBC = faces[0];
        PetscInt idRLBC = iT * numCells + faceLBC.cellIDR;
        addValuesBlock(jacobianAllGlobal, idRLBC, idRLBC, faceLBC.jacobianFluxRNeg[iT]);

        // --- 1.3 Flux Jacobian from the right boundary faces ---      
        const auto &faceRBC = faces[numCells];
        PetscInt idLRBC = iT * numCells + faceRBC.cellIDL;
        addValuesBlock(jacobianAllGlobal, idLRBC, idLRBC, faceRBC.jacobianFluxL[iT]);
    
        // --- 2.  Mass + chemical + joule source Jacaobian ---
        for (label i = 0; i < numCells; ++i)
        {
            PetscInt idC = iT * numCells + i;
            PetscInt idR = idC + 1;
            PetscInt idL = idC - 1;

            // --- 2.1 Mass + chemical source ---
            auto sumDiag = cells[i].massDiagnal + cells[i].jacobianCs[iT];
            addValuesBlock(jacobianAllGlobal, idC, idC, sumDiag);

            // --- 2.2 Joule source contribution ---
            if (i > 0) addValuesBlock(jacobianAllGlobal, idC, idL, jacobianFluxJouleLCell[iT][i]);
            if (i < numCells - 1) addValuesBlock(jacobianAllGlobal, idC, idR, jacobianFluxJouleRCell[iT][i]);
            addValuesBlock(jacobianAllGlobal, idC, idC, jacobianFluxJouleCCell[iT][i]);
        }
        // --- 3.  Harmonic balance source Jacobian ---
        for (label jT = 0; jT < numT; ++jT)
        {
            scalar Eij = harmMat(iT, jT);
            for (label i = 0; i < numCells; ++i)
            {
                PetscInt row = iT  * numCells + i;
                PetscInt col = jT  * numCells + i;
                jacobianHBs  = Eij * identity4;
                addValuesBlock(jacobianAllGlobal, row, col, jacobianHBs);
            }
        }
    }
    MatAssemblyBegin(jacobianAllGlobal, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacobianAllGlobal, MAT_FINAL_ASSEMBLY);
}


void FVM::Solver::getSlope()
{
    using namespace Const;
    using Tools::abs;
    using Tools::sign;

    // MUSCL-type slope limiter for computing a weighted average of left and right gradients
    auto slopeLimiter = [&](scalar fL, scalar fC, scalar fR, scalar dxL, scalar dxR) 
    {
        scalar sL = (fC - fL) / dxL;
        scalar sR = (fR - fC) / dxR;
        return (sign(sL) + sign(sR)) * abs(sL) * abs(sR) / (abs(sL) + abs(sR) + 1e-50);
    };

    for (label i = 0; i < numCells; ++i)
    {
        if (Const::isFirstOrder == "yes")
            break;

        Cell *cellL = &cells[i], *cellR = &cells[i + 1], *cellC = &cells[i];

        if (i == 0)
        {
            cellL = &cells[numCells + 1];
        }
        else
        {
            cellL = &cells[i - 1];
        }

        for (label iT = 0; iT < numT; ++iT)
        {
            cellC->sNe(iT)   = slopeLimiter(cellL->Ne(iT), cellC->Ne(iT), cellR->Ne(iT), faces[i].dist, faces[i + 1].dist);
            cellC->sNi(iT)   = slopeLimiter(cellL->Ni(iT), cellC->Ni(iT), cellR->Ni(iT), faces[i].dist, faces[i + 1].dist);
            cellC->sEe(iT)   = slopeLimiter(cellL->Ee(iT), cellC->Ee(iT), cellR->Ee(iT), faces[i].dist, faces[i+1].dist);
            cellC->sPhi(iT)  = slopeLimiter(cellL->Phi(iT), cellC->Phi(iT), cellR->Phi(iT), faces[i].dist, faces[i+1].dist);
        }
    }
    // Set the slopes at the boundary cells
    for (label iT = 0; iT < numT; ++iT)
    {
        cells[numCells].sNe(iT)      = 0.0;
        cells[numCells].sNi(iT)      = 0.0;
        cells[numCells].sEe(iT)      = 0.0;
        cells[numCells].sPhi(iT)     = 0.0;

        cells[numCells + 1].sNe(iT)  = 0.0;
        cells[numCells + 1].sNi(iT)  = 0.0;
        cells[numCells + 1].sEe(iT)  = 0.0;
        cells[numCells + 1].sPhi(iT) = 0.0;
    }
}

void FVM::Solver::evolve()
{
    using Const::numCells;
    for (auto &face : faces)
    {
        // Obtain the values at right and left sides of the faces
        face.interpolate();
        if (face.faceID == 0)
        {
            face.getFluxLeftBC();
        }
        else if (face.faceID == numCells)
        {
            face.getFluxRightBC();  
        }
        else
        {
            face.getFlux();
        }  
    }
}

void FVM::Solver::updateFluidExplicit(const label iRK)
{
    using namespace Const;
    using namespace Tools;

    for (label i = 0; i < numCells; ++i)
    {
        const scalar dtOverVol = cells[i].dt / cells[i].vol;
        for (label iT = 0; iT < numT; ++iT)
        {
            const scalar TeV       = KToeV(EeToTe(cells[i].Ne(iT) * nRef, cells[i].Ee(iT) * EeRef));
            cells[i].kl(iT)        = chemSets.interpolateChem(TeV) / klRef;


            cells[i].cS(iT)        = cells[i].kl(iT) * cells[i].Ne(iT) * N;
            cells[i].Je(iT)        = (faces[i].fluxNe(iT) + faces[i + 1].fluxNe(iT)) / 2.0;
            cells[i].Ji(iT)        = (faces[i].fluxNi(iT) + faces[i + 1].fluxNi(iT)) / 2.0;
           
            cells[i].ResFluxNe(iT) = alpha[iRK] * (faces[i].fluxNe(iT) - faces[i + 1].fluxNe(iT)) * dtOverVol;
            cells[i].ResFluxNi(iT) = alpha[iRK] * (faces[i].fluxNi(iT) - faces[i + 1].fluxNi(iT)) * dtOverVol;
            cells[i].ResFluxEe(iT) = alpha[iRK] * (faces[i].fluxEe(iT) - faces[i + 1].fluxEe(iT)) * dtOverVol;
            cells[i].ResFluxEe(iT) += alpha[iRK] * (faces[i].fluxEeJoule(iT) - faces[i + 1].fluxEeJoule(iT))
                                   * cells[i].Phi(iT) * dtOverVol;

            cells[i].ResFluxNe(iT) += alpha[iRK] * cells[i].kl(iT) * cells[i].Ne(iT) * cells[i].dt * N;
            cells[i].ResFluxNi(iT) += alpha[iRK] * cells[i].kl(iT) * cells[i].Ne(iT) * cells[i].dt * N;
            cells[i].ResFluxEe(iT) -= alpha[iRK] * cells[i].kl(iT) * cells[i].Ne(iT) * cells[i].dt * N * Hl;

        }

        cells[i].ResFluxNe -= alpha[iRK] * harmMat * cells[i].Ne * cells[i].dt;
        cells[i].ResFluxNi -= alpha[iRK] * harmMat * cells[i].Ni * cells[i].dt;
        cells[i].ResFluxEe -= alpha[iRK] * harmMat * cells[i].Ee * cells[i].dt;
        
        cells[i].Ne        += cells[i].ResFluxNe;
        cells[i].Ni        += cells[i].ResFluxNi;
        cells[i].Ee        += cells[i].ResFluxEe;
      
        for (label iT = 0; iT < numT; ++iT)
        {
            cells[i].Te(iT) = EeToTeND(cells[i].Ne(iT), cells[i].Ee(iT));
        }
    }

}


void FVM::Solver::updateFluidHBFCI()
{
    using namespace Const;
    using namespace Tools;

    // Create solver
    KSP HBFCISolver;
    // Create pre-conditioner
    PC HBFCIPC;
    KSPCreate(PETSC_COMM_SELF, &HBFCISolver); 

    // Set coefficient matrix and pre-conditioner matrix                      
    KSPSetOperators(HBFCISolver, jacobianAllGlobal, jacobianAllGlobal); 
    
    // ---------- LU decomposition ------------
    // KSPSetType(HBFCISolver, KSPPREONLY);
    // KSPGetPC(HBFCISolver, &HBFCIPC);
    // PCSetType(HBFCIPC, PCLU);

    // ---------------GMRES -------------------
    KSPSetType(HBFCISolver, KSPGMRES);                                 
    KSPGetPC(HBFCISolver, &HBFCIPC);
    PCSetType(HBFCIPC, PCHYPRE);
    PCHYPRESetType(HBFCIPC,"boomeramg");   
    KSPSetTolerances(HBFCISolver,
                 1e-8,   // rtol
                 1e-10,   // abstol
                 PETSC_DEFAULT,     // dtol
                 500);   // max_iter

    KSPSetFromOptions(HBFCISolver);                                
    // Solving: x = A⁻¹·b
    KSPSolve(HBFCISolver, rhsGlobal, incrementGlobal); 

    // Update the solution
    for (label iT = 0; iT < numT; ++iT)
    {
        for (label i = 0; i < numCells; ++i)
        {   PetscInt indices[4];
            indices[0] = 4 * iT * numCells + 4 * i + 0; 
            indices[1] = 4 * iT * numCells + 4 * i + 1;
            indices[2] = 4 * iT * numCells + 4 * i + 2;
            indices[3] = 4 * iT * numCells + 4 * i + 3;

            PetscScalar incrementW[4];

            VecGetValues(incrementGlobal, 4, indices, incrementW);

            cells[i].Ne(iT)  = cells[i].NeOld(iT)  + eps * incrementW[0];
            cells[i].Ni(iT)  = cells[i].NiOld(iT)  + eps * incrementW[1];
            cells[i].Ee(iT)  = cells[i].EeOld(iT)  + eps * incrementW[2];
            cells[i].Phi(iT) = cells[i].PhiOld(iT) + eps * incrementW[3];

            // Update the electron temperature
            cells[i].Te(iT)  = EeToTeND(cells[i].Ne(iT), cells[i].Ee(iT));
        }

        // Update the electric field
        const auto &faceLBC = faces[0];
        const auto &faceRBC = faces[numCells];
        for (label i = 0; i < numCells; ++i)
        {
            if (i == 0)
            {
                cells[i].Ec(iT) = - (cells[i].Phi(iT) - faceLBC.PhiL(iT)) / faceLBC.dist;
            }
            else if (i == numCells - 1)
            {
                cells[i].Ec(iT) = - (faceRBC.PhiR(iT) - cells[i].Phi(iT)) / faceRBC.dist;
            }
            else
            {
                cells[i].Ec(iT) = - (cells[i + 1].Phi(iT) - cells[i - 1].Phi(iT))
                                / (faces[i].dist + faces[i + 1].dist);
            }
        }
    }
    KSPDestroy(&HBFCISolver);  
}


void FVM::Solver::setBoundaryConditions(const Eigen::VectorXd &harmTime)
{
    using namespace Const;
    Cell *cellLBC = &cells[numCells + 1];
    Cell *cellRBC = &cells[numCells];
    
    for(label iT = 0; iT < numT; ++iT)
    {
        cellLBC->Ne(iT)  = cells[0].Ne(iT);
        cellLBC->Ni(iT)  = cells[0].Ni(iT);
        cellLBC->Ee(iT)  = cells[0].Ee(iT);
        cellLBC->Te(iT)  = cells[0].Te(iT);
        if (analysisMode == AnalysisMode::STEADY)
        {
            cellLBC->Phi(iT) = PhiRF;
        }
        else
        {
            cellLBC->Phi(iT) = PhiRF * std::sin(2.0 * pi * fRF * harmTime(iT) + phase); 
        }
        cellLBC->Ec(iT)  = cells[0].Ec(iT);
        

        cellRBC->Ne(iT)  = cells[numCells - 1].Ne(iT);
        cellRBC->Ni(iT)  = cells[numCells - 1].Ni(iT);
        cellRBC->Ee(iT)  = cells[numCells - 1].Ee(iT);
        cellRBC->Te(iT)  = cells[numCells - 1].Te(iT);
        cellRBC->Phi(iT) = 0.0;
        cellRBC->Ec(iT)  = cells[numCells - 1].Ec(iT);
    }
}

void FVM::Solver::sumForAve()
{
    using namespace Const;

    for (label i = 0; i < numCells; ++i)
    {
        const scalar dtWeight = cells[i].dt * tRef;
        cells[i].sumNe  += cells[i].Ne(0)  * nRef   * dtWeight;
        cells[i].sumNi  += cells[i].Ni(0)  * nRef   * dtWeight;
        cells[i].sumEe  += cells[i].Ee(0)  * EeRef  * dtWeight;
        cells[i].sumTe  += cells[i].Te(0)  * TeRef  * dtWeight;
        cells[i].sumPhi += cells[i].Phi(0) * phiRef * dtWeight;
        cells[i].sumE   += cells[i].Ec(0)  * ERef   * dtWeight;
    }
}

void FVM::Solver::updateAve()
{
    using namespace Const;
    const scalar fRFweight = fRF * fRef;
    for (label i = 0; i < numCells; ++i)
    {
        // Calculate the averaged values within a period
        cells[i].NeA  = cells[i].sumNe  * fRFweight;
        cells[i].NiA  = cells[i].sumNi  * fRFweight;
        cells[i].EeA  = cells[i].sumEe  * fRFweight;   
        cells[i].TeA  = cells[i].sumTe  * fRFweight;
        cells[i].PhiA = cells[i].sumPhi * fRFweight;
        cells[i].EA   = cells[i].sumE   * fRFweight;

        // Reset the summation to zero
        cells[i].sumNe  = 0.0;
        cells[i].sumNi  = 0.0;
        cells[i].sumEe  = 0.0; 
        cells[i].sumTe  = 0.0;  
        cells[i].sumPhi = 0.0; 
        cells[i].sumE   = 0.0;   

    }
}

void FVM::Solver::handleAveraging(bool& isAveraging, bool& isAveWritten)
{
    using namespace Const;
    if (!isAveraging && runCyc >= stopPeriod - 1.0)
    {
        std::cout << "Start averaging in last cycle at running cycle = " << runCyc << std::endl;
        isAveraging = true;
    }

    if (isAveraging)
    {
        sumForAve();
    }

    if (isAveraging && !isAveWritten && runCyc >= stopPeriod)
    {
        updateAve();
        writeFinalAverage();
        isAveWritten = true;

        std::cout << "Final average written at running cycle = " << runCyc << std::endl;
    }

}

void FVM::Solver::writeFourierCoefficientsHB()
{
    using namespace Const;

    if (step % writeStep != 0)  return;

    const std::string filename = outputDir + "/FourierCoefficients.dat";
    std::ofstream outFourier(filename, std::ios::app); 
    if (!outFourier) {
        std::cerr << "Cannot open " << filename << '\n';
        return;
    }

    // Write data for the current step into a "zone" section with the proper format
    outFourier << "ZONE T=\"Step_" << step << "\", I=" << numCells << ", F=POINT\n";

    for (label i = 0; i < numCells; i++) {
        const Eigen::VectorXd HatNe  = dftMat * cells[i].Ne * nRef;
        const Eigen::VectorXd HatNi  = dftMat * cells[i].Ni * nRef;
        const Eigen::VectorXd HatEe  = dftMat * cells[i].Ee * EeRef;
        const Eigen::VectorXd HatTe  = dftMat * cells[i].Te * TeRef;
        const Eigen::VectorXd HatPhi = dftMat * cells[i].Phi * phiRef;
        const scalar               x = (coords[i + 1].x + coords[i].x) * LRef * 0.5;

        outFourier << std::setprecision(10) << cells[i].cellID
                   << '\t' << x
                   << '\t'<< HatNe(0) 
                   << '\t' << HatNi(0)
                   << '\t' << HatEe(0)
                   << '\t' << HatTe(0) 
                   << '\t' << HatPhi(0);

        for (label nh = 1; nh <= numH; ++nh) 
        {
            outFourier << '\t' << HatNe(2 * nh)  << '\t' << HatNe(2 * nh - 1)
                       << '\t' << HatNi(2 * nh)  << '\t' << HatNi(2 * nh - 1)
                       << '\t' << HatEe(2 * nh)  << '\t' << HatEe(2 * nh - 1)
                       << '\t' << HatTe(2 * nh)  << '\t' << HatTe(2 * nh - 1)
                       << '\t' << HatPhi(2 * nh) << '\t' << HatPhi(2 * nh - 1);
        }

        outFourier << std::endl;
    }

    outFourier.close();
}


void FVM::Solver::writeUnsteadyFlowFieldHB()
{
    using namespace Const;
        
    if (step % writeStep != 0)  return;

    const std::string filename = outputDir + "/UnsteadyFlowFieldHB.dat";
    std::ofstream outSolution(filename, std::ios::app); 
    if (!outSolution) {
        std::cerr << "Cannot open " << filename << '\n';
        return;
    }

    Eigen::VectorXd Ne(Eigen::VectorXd::Zero(numCells));
    Eigen::VectorXd Ni(Eigen::VectorXd::Zero(numCells));
    Eigen::VectorXd Ee(Eigen::VectorXd::Zero(numCells));
    Eigen::VectorXd Te(Eigen::VectorXd::Zero(numCells));
    Eigen::VectorXd Phi(Eigen::VectorXd::Zero(numCells));
    Eigen::VectorXd cosPhase(Eigen::VectorXd::Zero(numH));
    Eigen::VectorXd sinPhase(Eigen::VectorXd::Zero(numH));

    for (label nI = 0; nI < nTI; ++nI)
    {
        scalar toverT = (1.0 * nI) / nTI; 
        scalar t = toverT * period * tRef;

        Eigen::VectorXd phaseT = angFreq * t;

        outSolution << "ZONE T=\"Step_" << step << "_t=" << toverT 
                    << "T\", I=" << numCells << ", F=POINT" << std::endl;


        for (label nh = 0; nh < numH; ++nh)
        {
            cosPhase(nh) = cos(phaseT(nh));
            sinPhase(nh) = sin(phaseT(nh));
        }

        for (label i = 0; i < numCells; ++i)
        {
            const Eigen::VectorXd HatNe  = dftMat * cells[i].Ne * nRef;
            const Eigen::VectorXd HatNi  = dftMat * cells[i].Ni * nRef;
            const Eigen::VectorXd HatEe  = dftMat * cells[i].Ee * EeRef;
            const Eigen::VectorXd HatTe  = dftMat * cells[i].Te * TeRef;
            const Eigen::VectorXd HatPhi = dftMat * cells[i].Phi * phiRef;

            Ne(i)  = HatNe(0);
            Ni(i)  = HatNi(0);
            Ee(i)  = HatEe(0);
            Te(i)  = HatTe(0);
            Phi(i) = HatPhi(0);

            for (label nh = 0; nh < numH; ++nh)
            {
                Ne(i)  += HatNe(2 * nh + 1) * sinPhase(nh) + HatNe(2 * nh + 2) * cosPhase(nh);
                Ni(i)  += HatNi(2 * nh + 1) * sinPhase(nh) + HatNi(2 * nh + 2) * cosPhase(nh);
                Ee(i)  += HatEe(2 * nh + 1) * sinPhase(nh) + HatEe(2 * nh + 2) * cosPhase(nh);
                Te(i)  += HatTe(2 * nh + 1) * sinPhase(nh) + HatTe(2 * nh + 2) * cosPhase(nh);
                Phi(i) += HatPhi(2 * nh + 1) * sinPhase(nh) + HatPhi(2 * nh + 2) * cosPhase(nh);
            }
            const scalar x = (coords[i + 1].x + coords[i].x) * LRef * 0.5;
            outSolution << std::setprecision(10)
                        << cells[i].cellID << '\t'
                        << x << '\t'
                        << Ne(i) << '\t'
                        << Ni(i) << '\t'
                        << Ee(i) << '\t'
                        << Te(i) << '\t'
                        << Tools::KToeV(Te(i)) << '\t'
                        << Phi(i) << std::endl;
        }
    }
    outSolution.close();
}



void FVM::Solver::writeIterSolution()
{
    using namespace Const;

    if (step % writeStep != 0)  return;
    
    const std::string filename = outputDir + "/iterSolution.dat";
    std::ofstream outFile(filename, std::ios::app);  
    if (!outFile) {
        std::cerr << "Cannot open " << filename << '\n';
        return;
    }


    for (label iT = 0; iT < Const::numT; iT++)
    {
        outFile << "ZONE T=\"Step_" << step << "_nT=" << iT << "\",I=" << numCells << ", F=POINT" << std::endl;

        for (label i = 0; i < numCells; i++)
        {
            const scalar x = (coords[i + 1].x + coords[i].x) * LRef * 0.5;
            outFile << std::setprecision(8)
                    << cells[i].cellID << '\t'
                    << x << '\t' 
                    << cells[i].Ne(iT) * nRef    << '\t' 
                    << cells[i].Ni(iT) * nRef    << '\t' 
                    << cells[i].Ee(iT) * EeRef   << '\t'
                    << cells[i].Te(iT) * TeRef   << '\t' 
                    << cells[i].Phi(iT) * phiRef << std::endl; 
        }
    }
    outFile.close();
}


void FVM::Solver::gridMetrics()
{
    using namespace Const;

    // Compute cell volumes
    for(label i = 0; i < numCells; ++i)
    {
        cells[i].vol = coords[i + 1].x - coords[i].x; 
    }
    // cell volume for the ghost cell(face element) at the right boundary
    cells[numCells].vol = 0.0; 

    // cell volume for the ghost cell(face element) at the left boundary
    cells[numCells  + 1].vol = 0.0;

    // Compute face distances
    for (label i = 0; i < numNodes; ++i) {
        if (i == 0) {
            faces[i].dist  = 0.5 * cells[i].vol;
            faces[i].distL = 0.0; 
            faces[i].distR = 0.5 * cells[i].vol;
        }
        else if (i == numNodes - 1) {
            faces[i].dist  = 0.5 * cells[i - 1].vol;
            faces[i].distL = 0.5 * cells[i - 1].vol;
            faces[i].distR = 0.0;
        }
        else {
            faces[i].dist  = 0.5 * (cells[i].vol + cells[i - 1].vol);
            faces[i].distL = 0.5 * cells[i - 1].vol;
            faces[i].distR = 0.5 * cells[i].vol;
        }
    }
}

void FVM::Solver::checkNegativeStates()
{
    using namespace Const;
    for (label i = 0; i < numCells; ++i)
    {
        for (label iT = 0; iT < numT; ++iT)
        {
            if (cells[i].Ne(iT) < 0.0 || 
                cells[i].Ni(iT) < 0.0 || 
                cells[i].Ee(iT) < 0.0 || 
                cells[i].Te(iT) < 0.0)
            {
                std::cerr << "Error: Negative state detected: "
                          << "cell=" << i 
                          << ", nT =" << iT
                          << ", Ne=" << cells[i].Ne(iT)
                          << ", Ni=" << cells[i].Ni(iT)
                          << ", Ee=" << cells[i].Ee(iT)
                          << ", Te=" << cells[i].Te(iT) 
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

void FVM::Solver::writeCellVolumePlot(const std::string& filename)
{
    using namespace Const;

    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "Error: Cannot open " << filename << " for writing cell volume plot.\n";
        return;
    }

    out << "TITLE = Cell Volume Distribution\n";
    out << "VARIABLES = CellID, Volume\n";

    // Zone 1: dimensional volume (with units)
    out << "ZONE T = \"Dimensional Volume\", F = POINT\n";
    for (label i = 0; i < numCells; ++i)
    {
        out << i << '\t' << std::setprecision(10) << cells[i].vol * LRef << '\n';
    }

    // Zone 2: nondimensional volume (original value)
    out << "\nZONE T = \"Nondimensional Volume\", F = POINT\n";
    for (label i = 0; i < numCells; ++i)
    {
        out << i << '\t' << std::setprecision(10) << cells[i].vol << '\n';
    }

    out.close();
    std::cout << "Cell volume plot written to: " << filename << "\n";
}




void FVM::Solver::readQuasiSteadySolution(const std::string& filePath)
{
    using namespace Const;
    std::ifstream inFile(filePath.c_str());

    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filePath << std::endl;
        return;
    }

    std::string line;

  
    while (std::getline(inFile, line)) {
        if (line.find("ZONE") != std::string::npos) {
            break;  
        }
    }

    label iT = 0; 
    while (true)
    {

        if (line.find("ZONE") != std::string::npos)
        {

            size_t pos = line.find("T=\"TimeInstant_");
            if (pos != std::string::npos)
            {
                std::string timeStr = line.substr(pos + 15); 
                iT = std::stoi(timeStr.substr(0, timeStr.find("\"")));  
            }


            if (!std::getline(inFile, line)) break;
        }


        std::istringstream iss(line);
        label cellID;
        double Ne_val, Ni_val, Te_val, Phi_val;

        if (!(iss >> cellID >> Ne_val >> Ni_val >> Te_val >> Phi_val)) {

            if (!std::getline(inFile, line)) break;
            continue;
        }


        for (label i = 0; i < numCells; i++)
        {
            if (cells[i].cellID == cellID)  
            {
                cells[i].Ne(iT)  = Ne_val;
                cells[i].Ni(iT)  = Ni_val;
                cells[i].Te(iT)  = Te_val;
                cells[i].Phi(iT) = Phi_val;
                break;
            }
        }

       
        if (!std::getline(inFile, line)) break;
    }

    inFile.close();
}



void FVM::Solver::writeResidual(scalar wallTime)
{
    using namespace Const;
    
    if (step % resWriteStep != 0)  return;

    const std::string filename = outputDir + "/ResidualHistory.dat";
    std::ofstream out(filename, std::ios::app); 
    if (!out) {
        std::cerr << "Cannot open " << filename << '\n';
        return;
    }

    out << step     << '\t'         
        << wallTime << '\t';
    if (analysisMode == AnalysisMode::DT)
    out << runCyc   << '\t';
    out << std::fixed << std::setprecision(8)
        << resNe(0) << '\t'
        << resNi(0) << '\t'
        << resEe(0) << '\t'
        << resTe(0) << '\n';
    out.close();

}

void FVM::Solver::initializeOutputFiles()
{
    using namespace Const;

    {
        std::ofstream out(outputDir + "/ResidualHistory.dat", std::ios::trunc);
        if (out) 
        {
            out << "TITLE = Residual History\n";
            out << "VARIABLES = Step, wallTime";
            if (analysisMode == AnalysisMode::DT)
                out << ", runCyc";
            out << ", log10(Res_Ne), log10(Res_Ni), log10(Res_Ee), log10(Res_Te)\n";
        }
        else 
        {
            std::cerr << "Failed to create ResidualHistory.dat\n";
        }
    }

    {
        std::ofstream out(outputDir + "/iterSolution.dat", std::ios::trunc);
        if (out) 
        {
            out << "TITLE= Solution Data\n";
            out << "VARIABLES= cellID, x, Ne, Ni, Ee, Te, Phi\n";
        } else 
        {
            std::cerr << "Failed to create iterSolution.dat\n";
        }
    }

    {
        if (analysisMode == AnalysisMode::DT 
            || (analysisMode == AnalysisMode::HB && numH == 0))
        {
            std::ofstream out(outputDir + "/NeNiEeTeAveTs.dat", std::ios::trunc);
            if (out) 
            {
                out << "variables = step, wallTime, runCyc, NeAve, NiAve, EeAve, TeAve\n";
            } else 
            {
                std::cerr << "Failed to create NeNiEeTeAveTs.dat\n";
            }
        }
    }

    {
        if (analysisMode == AnalysisMode::HB){
            std::ofstream out(outputDir + "/FourierCoefficients.dat", std::ios::trunc);
            if (out) 
            {
                out << "variables = cellID, x, Ne_mean, Ni_mean, Ee_mean, Te_mean, Phi_mean";
                for (label nh = 1; nh <= numH; ++nh) 
                {
                    out << ", Ne_real_" << nh << ", Ne_imag_" << nh
                        << ", Ni_real_" << nh << ", Ni_imag_" << nh
                        << ", Ee_real_" << nh << ", Ee_imag_" << nh
                        << ", Te_real_" << nh << ", Te_imag_" << nh
                        << ", Phi_real_" << nh << ", Phi_imag_" << nh;
                }
                out << '\n';
            } 
            else 
            {
                std::cerr << "Failed to create NeNiEeTeAveTs.dat\n";
            }
        }
    }

    {
        if (analysisMode == AnalysisMode::HB){
            std::ofstream out(outputDir + "/UnsteadyFlowFieldHB.dat", std::ios::trunc);
            if (out) 
            {
                out << "variables = cellID, x, nE, nI, Ee(J), Te(K), Te(eV), Phi\n";
            } 
            else 
            {
                std::cerr << "Failed to create NeNiEeTeAveTs.dat\n";
            }
        }
    }



}

void FVM::Solver::writeAverageNeNiTeTs(scalar wallTime)
{
    using namespace Const;

    if (step % writeStep != 0)  return;
    
    const std::string filename = outputDir + "/NeNiEeTeAveTs.dat";
    std::ofstream out(filename, std::ios::app); 
    if (!out) {
        std::cerr << "Cannot open " << filename << '\n';
        return;
    }


    const scalar invTotalCount = 1.0 / numCells;

    scalar sumNe = 0.0, sumNi = 0.0, sumEe = 0.0, sumTe = 0.0;


    for (label i = 0; i < numCells; ++i)
    {
        sumNe += cells[i].Ne(0);
        sumNi += cells[i].Ni(0);
        sumEe += cells[i].Ee(0);
        sumTe += cells[i].Te(0);
    }


    const scalar avgNe = sumNe * invTotalCount * nRef;
    const scalar avgNi = sumNi * invTotalCount * nRef;
    const scalar avgEe = sumEe * invTotalCount * EeRef;
    const scalar avgTe = sumTe * invTotalCount * TeRef;

 
    out << step << '\t'
        << std::fixed << std::setprecision(8) 
        << wallTime << '\t'
        << runCyc << '\t' 
        << std::scientific << std::setprecision(8)
        << avgNe << '\t'
        << avgNi << '\t'
        << avgEe << '\t'
        << avgTe << '\n';

    out.close();
}

void FVM::Solver::writeFinalAverage()
{
    using namespace Const;
    std::stringstream fileName;

    fileName << outputDir << "/primFinalAve_space.dat" << std::endl;

    std::string s;
    fileName >> s;

    std::ofstream out(s.c_str());


    out << "variables = X, NeA, NiA, EeA ,TeA(J), TeA(eV), PhiA, EA"
        << std::endl;

    for (label i = 0; i < numCells ; i++)
    {
        const scalar x = (coords[i + 1].x + coords[i].x) * LRef * 0.5;
        out << std::fixed << std::setprecision(10)
            << x << '\t'
            << cells[i].NeA << '\t'
            << cells[i].NiA << '\t'
            << cells[i].EeA << '\t'
            << cells[i].TeA << '\t'
            << Tools::KToeV(cells[i].TeA) << '\t'
            << cells[i].PhiA << '\t'
            << cells[i].EA << '\t'
            << std::endl;
    }

    out.close();

    std::cout << "Final average data written to: " << s << std::endl;

}


void FVM::Solver::calRes()
{
    using namespace Const;

    resNe.setZero();
    resNi.setZero();
    resEe.setZero();
    resTe.setZero();
    
    for (label i = 0; i < numCells; ++i)
    {
        for (label iT = 0; iT < numT; ++iT)
        {
            scalar neNew   = cells[i].Ne(iT);
            scalar neOld   = cells[i].NeOld(iT);
            scalar niNew   = cells[i].Ni(iT);
            scalar niOld   = cells[i].NiOld(iT);
            scalar eeNew   = cells[i].Ee(iT);
            scalar eeOld   = cells[i].EeOld(iT);
            scalar teNew   = cells[i].Te(iT);
            scalar teOld   = cells[i].TeOld(iT);

            scalar denomne = std::max(std::abs(neOld), 1e-10);  
            scalar denomni = std::max(std::abs(niOld), 1e-10);
            scalar denomee = std::max(std::abs(eeOld), 1e-10);  
            scalar denomte = std::max(std::abs(teOld), 1e-10);  

            resNe(iT)      += pow((neNew - neOld) / denomne, 2);
            resNi(iT)      += pow((niNew - niOld) / denomni, 2);
            resEe(iT)      += pow((eeNew - eeOld) / denomee, 2); 
            resTe(iT)      += pow((teNew - teOld) / denomte, 2);  

        }

    }
    for (label iT = 0; iT < numT; ++iT)
    {
        resNe(iT) = log10(sqrt(resNe(iT) / numCells) + 1e-20);
        resNi(iT) = log10(sqrt(resNi(iT) / numCells) + 1e-20);
        resEe(iT) = log10(sqrt(resEe(iT) / numCells) + 1e-20);
        resTe(iT) = log10(sqrt(resTe(iT) / numCells) + 1e-20);
    }
}

void FVM::Solver::infoRes()
{
    using namespace Const;
    

    if (step % printStep != 0)  return;



    calRes();

    std::cout << "step = " << step << '\t'          
              << std::scientific << std::setprecision(4);

    if (analysisMode == AnalysisMode::DT || (analysisMode == AnalysisMode::HB && numH == 0))       
        std::cout << "period = " << runCyc << '\t';

    std::cout << "Max Res(log10):  "
              << "Ne = " << resNe(0) << '\t'
              << "Ni = " << resNi(0) << '\t'
              << "Ee = " << resEe(0) << '\t'
              << "Te = " << resTe(0)
              << std::endl;
}

bool FVM::Solver::isExplicitDT() const
{
    using namespace Const;
    return analysisMode == AnalysisMode::DT 
        || (analysisMode == AnalysisMode::HB && numH == 0);
}




void FVM::Solver::iterateExplicit()
{
    using namespace Const;
    bool isAveraging         = false;     
    bool isAveWritten        = false;
    bool isWritingStepCal    = false;
    bool isResWritingStepCal = false; 
    bool isPrintStepCal      = false;   
    auto stepStart = std::chrono::high_resolution_clock::now();

    
    std::ofstream stopFile(outputDir +"/stop.dat", std::ios::trunc);
    if (!stopFile.is_open()) {
        std::cerr << "Error: Cannot create stop.dat file!" << std::endl;
        return;
    }
    stopFile << "0" << std::endl; 
    stopFile.close(); 

    std::cout << "Starting iterating..." << std::endl;
    while (
        stopMode == StopMode::TIME   ? runTime < stopTime :
        stopMode == StopMode::STEP   ? step    < stopStep :
        stopMode == StopMode::CYC ? runCyc  < stopPeriod :
        false
    )
    {
        std::ifstream checkStopFile(outputDir +"/stop.dat");
        if (checkStopFile.is_open()) 
        {
            label stopSignal = 0;
            checkStopFile >> stopSignal;
            checkStopFile.close();

            if (stopSignal == 1) { 
                std::cout << "Stop signal received, exiting loop." << std::endl;
                break; 
            }
        } 
        else 
        {
            std::cerr << "Error: Cannot open stop.dat for reading!" << std::endl;
        }

        // Get the time step for each cell
        getDt(isWritingStepCal, isResWritingStepCal, isPrintStepCal);

        // Store the conservative variables from the last step 
        initRK();

        // Runge-Kutta iteration
        label iRK = 0;
        while ( iRK < nRK )
        {
            // Get the slopes based on MUSCL limiter
            getSlope();

            // Get the fluxes of each face
            evolve();
            
            // Update electric potential and electric field in cells
            updatePhiFVM();
            
            // Update the flow field
            updateFluidExplicit(iRK);

            iRK++;
        }
        auto stepEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<scalar> wallTime = stepEnd - stepStart;

        step++;

        if (isExplicitDT())
        {
            harmTime[0] += dtminGlobal;
            runTime     += dtminGlobal;
            runCyc       = runTime / period;
        }

        if (isExplicitDT())
        {
            handleAveraging(isAveraging, isAveWritten);
        }



        // Output spaced-averaged values 
        if (isExplicitDT())
        {
            writeAverageNeNiTeTs(wallTime.count()); 
        }
  

        // Setup the boundary conditions
        setBoundaryConditions(harmTime);

        infoRes();

        writeResidual(wallTime.count());

        checkNegativeStates();

        writeIterSolution();


        if(analysisMode == AnalysisMode::HB)
        {
            writeFourierCoefficientsHB();
            writeUnsteadyFlowFieldHB();    

        }

    }
}


void initDebugLogs()
{
    if (!freopen("debug_log", "w", stdout)) {
        std::cerr << "Failed to redirect stdout to debug_log\n";
    }
    if (!freopen("debug_error.log", "w", stderr)) {
        std::cerr << "Failed to redirect stderr to debug_error.log\n";
    }
    std::cout.setf(std::ios::unitbuf);
}



void FVM::Solver::iterateHBFCI()
{
    using namespace Const;  
    const auto stepStart = std::chrono::high_resolution_clock::now();// 少用auto
 
    std::ofstream stopFile(outputDir +"/stop.dat", std::ios::trunc);
    if (!stopFile.is_open()) 
    {
        std::cerr << "Error: Cannot create stop.dat file!" << std::endl;
        return;
    }
    stopFile << "0" << std::endl; 
    stopFile.close(); 

    std::cout << "Starting iterating..." << std::endl;
    while (step < stopStep)
    {
        std::ifstream checkStopFile(outputDir +"/stop.dat");
        if (checkStopFile.is_open()) 
        {
            label stopSignal = 0;
            checkStopFile >> stopSignal;
            checkStopFile.close();

            if (stopSignal == 1) { 
                std::cout << "Stop signal received, exiting loop." << std::endl;
                break; 
            }
        } 
        else 
        {
            std::cerr << "Error: Cannot open stop.dat for reading!" << std::endl;
        }

        //- Get the time step for each cell
        std::cout << "Getting time step for each cell..." << std::endl;
        getDtau();

        //- Get the mass stiffness diagnal：vol / dtau
        std::cout << "Getting mass stiffness diagnal..." << std::endl;
        getMassDiagnal();

        //- Store the conservative variables from the last step 
        std::cout << "Storing conservative variables from the last step..." << std::endl;
        initRK();

        //- Runge-Kutta iteration
        label iRK = 0;
        while ( iRK < nRK )
        {
            //- Get the slopes based on MUSCL limiter
            std::cout << "Getting slopes for Runge-Kutta step " << iRK << "..." << std::endl;
            getSlope();

            //- Get the fluxes of each face
            std::cout << "Evolving fluxes for Runge-Kutta step " << iRK << "..." << std::endl;
            evolve();

            //- Get the global RHS vector
            std::cout << "Assembling global RHS vector for Runge-Kutta step " << iRK << "..." << std::endl;
            assembleGlobalVecRHS(iRK);

            //- Get the local flux Jacobian
            std::cout << "Assembling local flux Jacobian for Runge-Kutta step " << iRK << "..." << std::endl;
            assembleLocalFluxJacobian();

            //- Get the flux Joule Jacobian
            std::cout << "Getting flux Joule Jacobian for Runge-Kutta step " << iRK << "..." << std::endl;
            getFluxJouleJacobian();

            //- Get the chemical source term Jacobian 
            std::cout << "Getting chemical source term Jacobian for Runge-Kutta step " << iRK << "..." << std::endl;
            getCsJacobian();

            //- Get the global Jacobian
            std::cout << "Assembling global Jacobian for Runge-Kutta step " << iRK << "..." << std::endl; 
            assembleGlobalJacobian();

            //- Solve the linear system and update the flow variables
            std::cout << "Solving linear system for Runge-Kutta step " << iRK << "..." << std::endl;
            updateFluidHBFCI();

            iRK++;
        }
        auto stepEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<scalar> wallTime = stepEnd - stepStart;

        step++;

        //- Setup the boundary conditions
        std::cout << "Setting up boundary conditions..." << std::endl;
        setBoundaryConditions(harmTime);

        //- Calculate the residuals and print them
        std::cout << "Calculating and printing residuals..." << std::endl;
        infoRes();

        //- Write the residuals to the file
        writeResidual(wallTime.count());

        //- Check if negative states in the flow variables exist
        checkNegativeStates();

        //- Write the Quasi-steady solution to the output file
        writeIterSolution();

        if(analysisMode == AnalysisMode::HB)
        {
            //- Write the Fourier coefficients to the output file
            writeFourierCoefficientsHB();

            //- Write the reconstructed unsteady flow field to the output file
            writeUnsteadyFlowFieldHB();    
        }
    }
}

void FVM::Solver::iterateHBPCI1()
{

}

int main(int argc, char *argv[])
{
    using namespace Const;

    initDebugLogs();

    PetscInitialize(&argc, &argv, NULL, NULL);


    std::string caseFile;  

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-c" && i + 1 < argc) {
            caseFile = argv[i + 1];
            i++;
        }
    }   
    if (caseFile.empty()) {
        std::cerr << "Error: Missing -c <config_file>" << std::endl;
        PetscFinalize();
        return 1;
    }

    std::cout << "Loading case file: " << caseFile << std::endl;
    // Read Input Configuration File
    try {
            readConfigFile(caseFile);
        }
    catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }


    std::cout << "Scaling to dimensionless variables" << std::endl;
    // Reading the mesh file
    readMeshFile(meshFile);
    
    std::cout << "Writing mesh file" << std::endl;
    // Writing the mesh file
    writeNodeLinePlot("mesh_nodes_line.dat");

    std::cout << "Non-dimensionalizing the parameters and variables" << std::endl;
    // Non-dimensionalize the Parameters and Variables
    scaleToDimensionless();
    
    std::cout << "Initializing the FVM Solver" << std::endl;
    // Initialize the Solver
    FVM::Solver ccp;

    std::cout << "Building output files" << std::endl;
    // Building the output file
    ccp.initializeOutputFiles();

    std::cout << "Calculating the grid metrics" << std::endl;
    // Calculating the metrics of the grid
    ccp.gridMetrics();

    std::cout << "Writing cell volume" << std::endl;
    // Writing cell wolume 
    ccp.writeCellVolumePlot("mesh_cell_volume.dat");

    // Initializing the varaibles from default values or from the file
    std::cout << "Initializing the variables" << std::endl;
    ccp.initlize(); ///////!!!!

    // Initializing the harmonic matrix
    std::cout << "Initializing the matrices and vectors related to harmonic balance method" << std::endl;
    ccp.initHarmonicMat();

    // Setting up the poisson coefficient matrix
    if(Const::implicitScheme == ImplicitScheme::NO ||
       Const::implicitScheme == ImplicitScheme::PCI1)
    {
        std::cout << "Setting up thr poisson coefficient matrix: A" << std::endl;
        ccp.setupPoisson();
    }
    // Setup the boundary conditions
    std::cout << "Seting up the boundary condition" << std::endl;
    ccp.setBoundaryConditions(ccp.harmTime);///////!!!!


    if(analysisMode == AnalysisMode::HB)
    {
        if(implicitScheme == ImplicitScheme::PCI1)
        {
            std::cout << "Running HB implicitly coupling governing equations 1 2 3!" << std::endl;
            ccp.iterateHBPCI1();
        }
        else if (implicitScheme == ImplicitScheme::FCI)
        {
            std::cout << "Running HB implicitly coupling all governing equations!" << std::endl;
            ccp.iterateHBFCI();
        }
        else if (implicitScheme == ImplicitScheme::NO)
        {   
            std::cout << "Running HB explicitly!" << std::endl;
            ccp.iterateExplicit();
        }
        
    }
    else if (analysisMode == AnalysisMode::DT)
    {
        std::cout << "Running DT!" << std::endl;
        ccp.iterateExplicit(); 
    }
    
    std::cout << "Computation finished!🍽️😄" << std::endl;

    return 0;
}
