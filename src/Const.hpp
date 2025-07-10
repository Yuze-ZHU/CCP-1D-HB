#include <valarray>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "def.hpp"

#include "../eigen-3.3.9/Eigen/Dense"
using label = long long int;
using scalar = double;
using tensor1 = std::valarray<scalar>;
using tensor2 = std::valarray<std::valarray<scalar>>;
using std::vector;

namespace Const
{

    // ------------------------------------------------- //
    // Description:                                      //
    // Save initial condition, boundary condition        //
    // and gas properties                                //
    //                                                   //
    // Case:                                             //
    // Lymberopoulos 1D                                  //
    // doi: 10.1063/1.352926                             //
    // ------------------------------------------------- //

    //-------------- Input File ------------------------

    // - Initialized file
    extern std::string initFile;

    // - Rate coefficient file
    extern std::string chemCSV;

    // - Mesh file
    extern std::string meshFile;


    //------------- Mesh Infomation -------------------
    // - Length, m
    extern scalar L;

    //- Coordinates of the nodes
    extern vector<NODE> coords;

    //- Number of nodes
    extern label numNodes;

    //- Number of cells
    extern label numCells;


    //------------- Calculation Setup -------------------

    //- Number of harmonics
    extern label numH;

    //- Number of time instants
    extern label numT;

    
    //- Blocksize
    extern label blockSizeS;
    extern label blockSizeM;
    extern label blockSizeL;

    //- Voltage amplitude of RF, V
    extern scalar PhiRF;

    //- Frequency of RF, Hz
    extern scalar fRF;

    // - Initial phase
    extern scalar phase;
    
    //- Period of RF, s
    extern scalar period;

    //- Runge-Kutta Time Step
    extern label nRK;

    //- Runge-Kutta coefficients
    extern std::vector<scalar> alpha;

    //- CFL number
    extern scalar CFL;
    
    //- Relaxation factor for updating solution implicitly
    extern scalar eps;
    
    //- Number of loops for Block Jacobi iteration
    extern label nBJ;
    
    //- Error tolerance for Block-Jacobi iteration
    extern scalar tolNe;
    extern scalar tolNi;
    extern scalar tolEe;

    //- Number of variables treated implicitly
    extern label numVar;

    //- Initial Condition
    extern scalar Ni0;
    extern scalar Ne0;
    extern scalar Te0;
    extern scalar Phi0;

    //- Nondimensionalization
    extern scalar LRef;    // m
    extern scalar TeRef;   // K 
    extern scalar tRef;   // s
    extern scalar fRef;   // Hz
    extern scalar EeRef;   // Jm^-3
    extern scalar nRef;    // m^-3
    extern scalar phiRef;  // V
    extern scalar DeRef;   // m^2/s
    extern scalar DiRef;   // m^2/s
    extern scalar muERef;  // m^2/s/V
    extern scalar muIRef;  // m^2/s/V
    extern scalar eRef;    // C
    extern scalar klRef;   // m^3/s
    extern scalar HlRef;   // J
    extern scalar epsRef;  // F m^-1
    extern scalar cSRef;
    extern scalar JeRef;
    extern scalar JiRef;
    extern scalar ERef;
    extern scalar omegaRef; // angular frequency of RF, rad/s

    //- Is first order or not
    extern std::string isFirstOrder;

    //- Option for the constranit of Je of the electon flux
    extern std::string isJeConstraint; // YES or NO

    //- Types of analysis modes:
    extern AnalysisMode analysisMode; //STEADY, DT, HB

    //- Options for initializing the solution:
    extern InitOption initOption; // SCRATCH, CONTINUE, RESTART, FROMSTEADY

    //- Implicit scheme for pseudo-time stepping:
    extern ImplicitScheme implicitScheme; // PCI1

    //- Kind of time step used in the simulation:
    extern TimeStepType timeStepType; 

    //- Conditions to stop the simulation:
    extern StopMode stopMode;

    //- Types of output-write mdoe
    extern WriteMode  writeMode;

    //- Types of residual write mode
    extern ResWriteMode resWriteMode;

    //- Types of print mdoe
    extern PrintMode  printMode;


    //-------------------- Macro Properity ------------------------

    //- PI
    extern scalar pi;

    // - Boltzmann constant, J/K
    extern scalar kB;

    // - Elementary charge, C
    extern scalar e;

    // - Vacuum permittivity, F/m
    extern scalar eps0;

    // - Electron voltage to Joule conversion factor
    extern scalar eV2J;

    // - Energy exchange per collision, J
    extern scalar Hl;

        // - Number density of Ar(background gas), m^-3
    extern scalar N;

    // - Electron diffusivity, m^2/s
    extern scalar De;

    // - Ion diffusivity, m^2/s
    extern scalar Di;

    // - Electron mobility, m^2/s/V
    extern scalar muE;

    // - Ion mobility, m^2/s/V
    extern scalar muI;

    // - Mass of electron, kg
    extern scalar me;


    //- Secondary electron emission coefficient
    extern scalar Gam;

    // - Pressure of background gas, torr
    extern scalar P;

    // - InitiaL number density of background gas, m^-3
    extern scalar N0;

    // - Electron diffusivity times number density of background gas
    extern scalar De0;

    // - Ion diffusivity times number density of background gas
    extern scalar Di0;

    // - Electron mobility times number density of background gas
    extern scalar muE0;

    // - Ion mobility times number density of background gas
    extern scalar muI0;

    


    //------------------ Stop and Postprocess ---------------------
    //- Stop Time
    extern scalar stopTime;

    //- Stop Step
    extern label stopStep;

    //- Stop Step W
    extern label stopPeriod;


    //- Write Interval
    extern scalar writeCyc;
    extern label  writeStep;

    //- Write Interval for residual
    extern scalar resWriteCyc;
    extern label  resWriteStep;
    
    //- Print Interval 
    extern scalar printCyc;
    extern label  printStep;

    //-Number of Cells for Time Series
    extern label nCellTs;
    extern std::vector<label> cellIDTs;

    //- Number of time instants used for reconstruction
    extern label nTI;

    void readConfigFile(const std::string &filename);
    void readMeshFile(const std::string &meshFile);
    void scaleToDimensionless();
    void writeNodeLinePlot(const std::string& filename);
    std::string toString(AnalysisMode mode);
    std::string toString(InitOption option);
    std::string toString(ImplicitScheme scheme);
    std::string toString(TimeStepType type);
    std::string toString(StopMode mode);
    std::string toString(WriteMode mode);
    std::string toString(ResWriteMode mode);
    std::string toString(PrintMode mode);
}