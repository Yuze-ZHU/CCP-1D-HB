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

#include <valarray>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include "def.hpp"
#include <petsc.h>
#include "../eigen-3.3.9/Eigen/Dense"
#include "../eigen-3.3.9/Eigen/Sparse"
// Typedefs

using label = long long int;
using scalar = double;
using tensor1 = std::valarray<scalar>;
using tensor2 = std::valarray<std::valarray<scalar>>;
using std::vector;
using RowMajorMatrixXd = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using DiagnalMatrixXd4 = Eigen::Matrix<scalar, 4, 4, Eigen::RowMajor>;
using DiagnalMatrixXd3 = Eigen::Matrix<scalar, 3, 3, Eigen::RowMajor>;


namespace Tools
{
    // Mathematic

    scalar abs(const scalar num);

    tensor1 sign(const tensor1 &array);
    
    scalar sign(const scalar &value);

    scalar max(const scalar a, const scalar b);

    scalar min(const scalar a, const scalar b);

    // Basic Kinetic Formula

    scalar KToeV(const scalar K);

    scalar eVtoK(const scalar eV);

    scalar TeToEe(const scalar Ne, const scalar Te);

    scalar TeToEeND(const scalar Ne, const scalar Te);

    scalar EeToTe(const scalar Ne, const scalar Ee);

    scalar EeToTeND(const scalar Ne, const scalar Ee);

    std::string toString(double value, int precision);

    std::string toStringFill0(double value);

    std::vector<PetscScalar> toRowMajor(const Eigen::MatrixXd& M);
}


namespace FVM
{

    // Basic Element for Finite Volume Element

    struct ChemRate
    {
        double energy;
        double cRate;
    };
    
    class ChemSet
    {

    public:
        // Member Variables

        //- Sets of Chemial Reaction Rate 
        vector<ChemRate> chemSets;

        // Constructor

        //- Constructed by .csv file
        explicit ChemSet(const std::string csvName);

        // Destructor
        ~ChemSet();

        // Member function
        scalar interpolateChem(const scalar Te);

    };


    class Cell
    {
    public:
        // Member Variables

        // Geometry
        // - Volume of the cell
        scalar vol;

        //- Cell ID
        label cellID;

        //- Physical time step
        scalar dt;

        //- Pseudo time step
        scalar dtau;

        //- Sum the values within the last cycle
        scalar sumNe, sumNi, sumEe, sumTe, sumPhi, sumE;

        //- Averaged values within the last cycle
        scalar NeA, NiA, EeA, TeA, PhiA, EA;

        //- Rate coefficient 
        Eigen::VectorXd kl;

        // Conservative variables
        Eigen::VectorXd Ne, Ni, Ee, Phi;

        // Primitive variables
        Eigen::VectorXd Te, Ec;

        //- Other variables
        Eigen::VectorXd Je, Ji, cS;

        //- Conservative variables at the previous time step
        Eigen::VectorXd NeOld, NiOld, EeOld, PhiOld;

        //- Primitive variables at the previous time step
        Eigen::VectorXd TeOld;

        //- Residual flux of conservative variables(including flux and source term, updating the conservative variables) 
        Eigen::VectorXd ResFluxNe, ResFluxNi, ResFluxEe, ResFluxPhi;
        
        //- Slope of conservative variables
        Eigen::VectorXd sNe, sNi, sEe, sPhi;

        //- Stiffness matrix
        RowMajorMatrixXd massDiagnal;

        //- Chemical source Jacobian
        std::vector<RowMajorMatrixXd> jacobianCs;
        
        // Constructor
        //- Construct by id
        explicit Cell(const label i);

        // Destructor
        ~Cell();
    };

    class Parcel
    {
    };

    class Face
    {

    public:
        // Member Variables

        //- Cells at left and right side
        const Cell &cellL, &cellR;

        //- Face ID
        label faceID;

        //- Cell ID
        label cellIDL, cellIDR;

        //- Distance between the centers of the left and right cells
        scalar dist;

        //- Distance between the center of the left cell and the face
        scalar distL;

        //- Distance between the center of the right cell and the face
        scalar distR;

        //- Conservative and primitive variables at left and right sides
        Eigen::VectorXd NeL, NeR;
        Eigen::VectorXd NiL, NiR;
        Eigen::VectorXd EeL, EeR;
        Eigen::VectorXd TeL, TeR;
        Eigen::VectorXd PhiL, PhiR;
        Eigen::VectorXd Ef;

        //- Flux for number density, electron energy and electric field
        Eigen::VectorXd fluxNe, fluxNi;
        Eigen::VectorXd fluxEe;
        Eigen::VectorXd fluxEeJoule;
        Eigen::VectorXd fluxPhi;

        //- Flux Jacobian
        std::vector<RowMajorMatrixXd> jacobianFluxL;
        std::vector<RowMajorMatrixXd> jacobianFluxR;
        std::vector<RowMajorMatrixXd> jacobianFluxLNeg;
        std::vector<RowMajorMatrixXd> jacobianFluxRNeg;

        Eigen::VectorXd jacobianFluxEeL;
        Eigen::VectorXd jacobianFluxEeR;

        // Constructor
        //- Construct by cells
        explicit Face(const Cell &Lcell, const Cell &Rcell);

        // Destructor
        ~Face();

        // Member function

        //- Interpolate to face
        void interpolate();

        //- Calculate diffusion term interface
        scalar getDiff(const scalar pLcell, const scalar pRcell);

        //- Calculate upwind term at cell interface for electron
        scalar getUpwind(const scalar pLcell, const scalar pRcell, const label iT);

        //- Calculate downwind term at cell interface for electron
        scalar getDownwind(const scalar pLcell, const scalar pRcell, const label iT);

        //- Calculate flux
        void getFlux();
        void getFluxLeftBC();
        void getFluxRightBC();

        //- Calculate flux Jacobian 
        void getFluxJacobianFCI();
        void getFluxJacobianLeftBCFCI();
        void getFluxJacobianRightBCFCI();

        void getFluxJacobianPCI1();
        void getFluxJacobianLeftBCPCI1();
        void getFluxJacobianRightBCPCI1();

        void getFluxJacobianPCI2();
        void getFluxJacobianLeftBCPCI2();
        void getFluxJacobianRightBCPCI2();

        void getFluxJacobianEe();
        void getFluxJacobianLeftBCEe();
        void getFluxJacobianRightBCEe();

        //- Get the negative value of flux jaocobian
        void getFluxJacobianNeg();
    };

    class Solver
    {
    public:
        // Member Variables

        //- Iteration step
        label step;

        //- Run time
        scalar runTime;
        
        //- Run period
        scalar runCyc;

        //- Time step
        scalar dtminGlobal;

        //- Is negative 
        bool isWarning;

        //- Output directory
        std::string outputDir;

        //- Chemical reacting set
        ChemSet chemSets;

        //- Basic FVM for cell
        vector<Cell> cells;

        //- Basic FVM for face
        vector<Face> faces;

        //- Residuals for output
        Eigen::VectorXd resNe, resNi, resEe, resTe;

        // Matrices  and vectors related to Harmonic Balance
        //- Time spectral source matrix, E
        RowMajorMatrixXd harmMat;

        //- DFT matrix, D
        RowMajorMatrixXd dftMat;

        //- Inverse of DFT matrix, D^-1
        RowMajorMatrixXd dftMatInv;
 
        //- Time instants
        Eigen::VectorXd harmTime;

        //- Angular frequencies
        Eigen::VectorXd angFreq;

        //- Harmonic balance source Jacobian
        RowMajorMatrixXd jacobianHBs;

        //- FluxJoule Jacobian(Joule source term)
        std::vector<std::vector<RowMajorMatrixXd>> jacobianFluxJouleLCell;
        std::vector<std::vector<RowMajorMatrixXd>> jacobianFluxJouleCCell;
        std::vector<std::vector<RowMajorMatrixXd>> jacobianFluxJouleRCell;        

        // Module for solving Poisson's equation using the FVM framework: A·x = b
        //- Coefficient matrix A for the discretized Poisson's equation
        Mat poissonMat;

        //- Solution vector, containing the electric potential (phi) values at each cell
        Vec phiVec;

        //- Right-hand side vector, including source terms and boundary conditions
        Vec phiRhsVec;

        //- Global Jacobian matrix including flux, chemical, and harmonic balance (HB) source terms across all time instants
        Mat jacobianAllGlobal;

        //- Global solution increment vector (ΔW), containing the implicitly-coupled updates 
        //  for all cells across all time instants in the harmonic balance system
        Vec incrementGlobal;

        //- Global right-hand side (RHS) vector of the coupled system, including flux, source, and HB terms
        Vec rhsGlobal;


        //- Global Jacobian matrix of fluxEe
        Mat jacobianEeGlobal;

        Vec incrementEe;

        //- Global right-hand side (RHS) vector of the PCI2 system, only including fluxEe;
        Vec rhsEe;

        //
        // Constructor
        Solver();

        // Destructor
        ~Solver(); 

        // Member function
        //- Add values to a block matrix
        inline void addValuesBlock(Mat &mat, PetscInt row, PetscInt col,
                            const RowMajorMatrixXd &block);

        //- Initialize the flow region
        void initlizeFluid();

        //- Initialize the Harmonic Matrix
        void initHarmonicMat();

        //- Calculate the metrics for the grid
        void gridMetrics();

        //- Read solution from the initial file
        void readQuasiSteadySolution(const std::string& filePath);

        //- Calculate time step
        void getDt(bool& isWritingStepCal, bool& isResWritingStepCal, bool& isPrintStepCal);
        void getDtau();
        //- Initialize Runge-Kutta
        void initRK();

        //- Calculate the slope of the variables
        void getSlope();

        //- Evolution
        void evolve();

        //- Solving
        void iterateExplicit();
        void iterateHBPCI1();
        void iterateHBPCI2();
        void iterateHBFCI();

        //- Update fluid field
        void updateFluidExplicit(const label iRK);
        void updateNeImplicitEK(const label iRK);
        void updateNiEeImplicitE(const label iRK);
        void updateFluidImplicitE(const label iRK);
        void updateFluidHBFCI();
        void updateFluidHBPCI1();
        void updateFluidHBPCI2();
        void updateEeHBPCI2(const label iRK);

        //- Get the Joule flux source Jacobian
        void getFluxJouleJacobianFCI();
        void getFluxJouleJacobianPCI1();
    
        //- Get chemical source Jacobian
        void getCsJacobianFCI();
        void getCsJacobianPCI1();
        void getCsJacobianPCI2();

        //- Get mass stiffness matrix
        void getMassDiagnal(); // used in PCI1, PCI2, FCI

        //- Assembel global RHS vector
        void assembleGlobalVecRHSFCI(const label iRK);
        void assembleGlobalVecRHSPCI1(const label iRK);
        void assembleGlobalVecRHSPCI2(const label iRK);
        void assembleGlobalVecRHSEe(const label iRK);

        //- Assemble global jacobian using FDM
        void assembleGlobalJacobianFDMForFirstCellFCI(const label iT, const label iC);
        void dumpLocalBlocksFromGlobalJ(const label iT, const label iC);

        //- Assembel local flux Jacobian matrix
        void assembleLocalFluxJacobianFCI();
        void assembleLocalFluxJacobianPCI1();
        void assembleLocalFluxJacobianPCI2();
        void assembleLocalFluxJacobianEe();


        //- Assemble global Jacobian matrix
        void assembleGlobalJacobian();
        void assembleGlobalJacobianPCI2();
        void assembleGlobalJacobianEe();

        //- Update Phi for electric potential
        void setupPoisson();
        void updatePhiFVM();

        //- Boundary condition setup
        void setBoundaryConditions();

        //- Sumation for Average
        void sumForAve();
        void updateAve();
        void handleAveraging(bool& isAveraging, bool& isAveWritten);

        //- Check
        void checkNegativeStates();

        //- Load init file
        void loadInitFile(const std::string &initFile);

        //- Write
        void initializeOutputFiles();
        void writeCellVolumePlot(const std::string& filename);
        void writeFourierCoefficientsHB();
        void writeUnsteadyFlowFieldHB();
        void writeIterSolution();
        void writeFinalSolution();
        void writeResidual(scalar wallTime);
        void writeAverageNeNiTeTs(scalar wallTime);
        void writeFinalAverage();

        //- Print residual information
        void infoRes();
        void calRes();

        // Calculation option
        bool isExplicitDT() const;



    };

}
