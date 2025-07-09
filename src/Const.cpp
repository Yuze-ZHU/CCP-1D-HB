#include <valarray>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Const.hpp"
#include <sstream>
#include <unordered_map>
#include "def.hpp"
#include <iomanip>

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
    std::string initFile;

    // - Rate coefficient file
    std::string chemCSV;

    // - Mesh file
    std::string meshFile;

    //------------- Mesh Infomation -------------------
    // - Length, m
    scalar L;

    //- Coordinates of the nodes
    vector<NODE> coords;

    //- Number of nodes
    label numNodes;

    //- Number of cells
    label numCells;

    //------------- Calculation Setup -------------------

    // - Number of harmonics
    label numH;

    // - Number of time instants
    label numT;

    //- Blocksize
    label blockSizeS;
    label blockSizeM;
    label blockSizeL;

    // - Voltage amplitude of RF, V
    scalar PhiRF;

    // - Frequency of RF, Hz
    scalar fRF;

    // - Initial phase
    scalar phase;
    
    // - Period of RF, s
    scalar period;


    // - Runge-Kutta Time Step
    label nRK;

    // - Runge-Kutta coefficients
    std::vector<scalar> alpha;

    // - CFL number
    scalar CFL;
    
    // - relaxation factor for implicit scheme
    scalar eps;
    
    // - Number of loops for Block Jacobi iteration
    label nBJ;
    
    // - Error tolerance for Block-Jacobi iteration
    scalar tolNe;
    scalar tolNi;
    scalar tolEe;

    // - Number of variables treated implicitly
    label numVar;

    // Initial Condition
    scalar Ni0;
    scalar Ne0;
    scalar Te0;
    scalar Phi0;

    // Nondimensionalization
    scalar LRef;    // m
    scalar TeRef;   // K 
    scalar tRef;   // s
    scalar fRef;   // Hz
    scalar EeRef;   // Jm^-3
    scalar nRef;    // m^-3
    scalar phiRef;  // V
    scalar DeRef;   // m^2/s
    scalar DiRef;   // m^2/s
    scalar muERef;  // m^2/s/V
    scalar muIRef;  // m^2/s/V
    scalar eRef;    // C
    scalar klRef;   // m^3/s
    scalar HlRef;   // J
    scalar epsRef;  // F m^-1 
    scalar cSRef;
    scalar JeRef;
    scalar JiRef;
    scalar ERef; // V/m
    scalar omegaRef; // angular frequency, rad/s

    // - Is first order or not
    std::string isFirstOrder;

    // - Option for the constranit of Je of the electon flux
    std::string isJeConstraint; // YES or NO

    // - Types of analysis modes:
    AnalysisMode analysisMode; //STEADY, DT, HB

    // - Options for initializing the solution:
    InitOption initOption; // SCRATCH, CONTINUE, RESTART, FROMSTEADY

    // - Implicit scheme for pseudo-time stepping:
    ImplicitScheme implicitScheme; // PCI1

    //- Kind of time step used in the simulation:
    TimeStepType timeStepType; 

    // - Conditions to stop the simulation:
    StopMode stopMode;
    
    //- Types of output-write mdoe
    WriteMode  writeMode;

    //- Types of residual write mdoe
    ResWriteMode resWriteMode;

    //- Types of print mdoe
    PrintMode  printMode;



    //-------------------- Macro Properity ------------------------

    //- PI
    scalar pi = M_PI;

    // - Boltzmann constant, J/K
    scalar kB;

    // - Elementary charge, C
    scalar e;

    // - Vacuum permittivity, F/m
    scalar eps0;

    // - Electron voltage to Joule conversion factor
    scalar eV2J;
    
    // - Energy exchange per collision, J
    scalar Hl;

    // - Number density of Ar(background gas), m^-3
    scalar N;

    // - Electron diffusivity, m^2/s
    scalar De;

    // - Ion diffusivity, m^2/s
    scalar Di;

    // - Electron mobility, m^2/s/V
    scalar muE;

    // - Ion mobility, m^2/s/V
    scalar muI;

    // - Mass of electron, kg
    scalar me;

    //- Secondary electron emission coefficient
    scalar Gam;

    // - Pressure of background gas, torr
    scalar P;

    // - InitiaL number density of background gas, m^-3
    scalar N0;

    // - Electron diffusivity times number density of background gas
    scalar De0;

    // - Ion diffusivity times number density of background gas
    scalar Di0;

    // - Electron mobility times number density of background gas
    scalar muE0;

    // - Ion mobility times number density of background gas
    scalar muI0;
    

    //------------------ Stop and Postprocess ---------------------
    //- Stop Time
    scalar stopTime;

    //- Stop Step
    label stopStep;

    //- Stop Step 
    label stopPeriod;

    //- Write Interval
    scalar writeCyc;
    label  writeStep;

    //- Write Interval for residual
    scalar resWriteCyc;
    label  resWriteStep;
    
    //- Print Interval 
    scalar printCyc;
    label  printStep;

    //-Number of Cells for Time Series
    label nCellTs;

    //- ID of Cell for the time series output
    std::vector<label> cellIDTs;

    // - Number of transient instants for output
    label nTI;

    void readConfigFile(const std::string &filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cout << "Error: Cannot open " << filename << std::endl;
            return;
        }

        std::string line;
        std::unordered_map<std::string, std::string> config;

        
        while (std::getline(file, line)) {
            
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);
                config[key] = value;
            }
        }

        auto trim = [](std::string &s) {
            s.erase(0, s.find_first_not_of(" \t\n\r")); 
            s.erase(s.find_last_not_of(" \t\n\r") + 1);  
        };

        if (config.count("initFile")) {
            initFile = config["initFile"];
            trim(initFile);}

        if (config.count("chemCSV")) {
            chemCSV = config["chemCSV"];
            trim(chemCSV);}

        if (config.count("meshFile")) {
            meshFile = config["meshFile"];
            trim(meshFile);}



        if (config.count("numH")) numH = std::stoi(config["numH"]);

        numT = numH * 2 + 1;

        if (config.count("PhiRF")) PhiRF = std::stod(config["PhiRF"]);
        
        if (config.count("fRF")) fRF = std::stod(config["fRF"]);
        
        period = 1.0 / fRF;

        if (config.count("phase")) phase = std::stod(config["phase"]);

        phase = phase / 180.0 * pi;

        if (config.count("nRK")) nRK = std::stoi(config["nRK"]);

        if (config.count("alpha")) {
            alpha.clear();
            std::stringstream ss(config["alpha"]);
            std::string val;
            while (std::getline(ss, val, ',')) {
                alpha.push_back(std::stod(val));
            }
        }

        if (config.count("CFL")) CFL = std::stod(config["CFL"]);

        if (config.count("eps")) eps = std::stod(config["eps"]);

        if (config.count("nBJ")) nBJ = std::stoi(config["nBJ"]);

        if (config.count("tolNe")) tolNe = std::stod(config["tolNe"]);

        if (config.count("tolNi")) tolNi = std::stod(config["tolNi"]);

        if (config.count("tolEe")) tolEe = std::stod(config["tolEe"]);

        if (config.count("Ni0")) Ni0 = std::stod(config["Ni0"]);

        if (config.count("Ne0")) Ne0 = std::stod(config["Ne0"]);

        if (config.count("Te0")) Te0 = std::stod(config["Te0"]);

        if (config.count("Phi0")) Phi0 = std::stod(config["Phi0"]);

        if (config.count("LRef")) LRef = std::stod(config["LRef"]);

        if (config.count("TeRef")) TeRef = std::stod(config["TeRef"]);

        if (config.count("nRef")) nRef = std::stod(config["nRef"]);

        if (config.count("phiRef")) phiRef = std::stod(config["phiRef"]);

        if (config.count("isFirstOrder")){
            isFirstOrder = config["isFirstOrder"];
            trim(isFirstOrder);}

        if (config.count("isJeConstraint")){
            isJeConstraint = config["isJeConstraint"];
            trim(isJeConstraint);}

        if (config.count("analysisMode")) {
            std::string modeStr = config["analysisMode"];
            trim(modeStr);

            if (modeStr == "STEADY") analysisMode = AnalysisMode::STEADY;
            else if (modeStr == "DT") analysisMode = AnalysisMode::DT;
            else if (modeStr == "HB") analysisMode = AnalysisMode::HB;
            else throw std::runtime_error("Invalid analysisMode: " + modeStr);
        }

        if (analysisMode != AnalysisMode::HB) 
        {
            if (numH != 0 || numT != 1)
            {
                std::cout << "[Config] analysisMode != HB: Overriding numH = " << numH 
                          << ", numT = " << numT << " => numH = 0, numT = 1\n";
                numH = 0;
                numT = 1;
            }
        }

        if (config.count("initOption")) {
            std::string initStr = config["initOption"];
            trim(initStr);

            if (initStr == "SCRATCH") initOption = InitOption::SCRATCH;
            else if (initStr == "CONTINUE") initOption = InitOption::CONTINUE;
            else if (initStr == "RESTART")  initOption = InitOption::RESTART;
            else if (initStr == "FROMSTEADY") initOption = InitOption::FROMSTEADY;
            else throw std::runtime_error("Unknown initOption: " + initStr);
        }

        if (config.count("implicitScheme")) {
            std::string schemeStr = config["implicitScheme"];
            trim(schemeStr);
            if (schemeStr == "NO") implicitScheme = ImplicitScheme::NO;
            else if (schemeStr == "PCI1") implicitScheme = ImplicitScheme::PCI1;
            else throw std::runtime_error("Unknown implicitScheme: " + schemeStr);
        }

        if(implicitScheme == ImplicitScheme::PCI1) {
            numVar = 3; 
        }

        if (config.count("timeStepType")) {
            std::string stepStr = config["timeStepType"];
            trim(stepStr);

            if (stepStr == "LOCAL") timeStepType = TimeStepType::LOCAL;
            else if (stepStr == "GLOBAL") timeStepType = TimeStepType::GLOBAL;
            else throw std::runtime_error("Unknown timeStepType: " + stepStr);
        } else {
            timeStepType = TimeStepType::GLOBAL; // 默认值
            std::cout << "timeStepType not set. Defaulting to GLOBAL.\n";
        }

        if (analysisMode == AnalysisMode::DT)
        {
            if (timeStepType == TimeStepType::LOCAL)
            {
                timeStepType = TimeStepType::GLOBAL;
                std::cout << "DT mode requires GLOBAL time stepping. Forcing timeStepType = GLOBAL.\n";
            }
            else
            {
                std::cout << "DT mode with GLOBAL time stepping confirmed.\n";
            }
        }
        else
        {
            std::cout << "Non-DT mode. timeStepType = "
                      << (timeStepType == TimeStepType::LOCAL ? "LOCAL" : "GLOBAL")
                      << " will be used (may have limited effect).\n";
        }

        

        if (config.count("stopMode")) {
            std::string stopStr = config["stopMode"];
            trim(stopStr);

            if (stopStr == "TIME")        stopMode = StopMode::TIME;
            else if (stopStr == "STEP")   stopMode = StopMode::STEP;
            else if (stopStr == "CYC")    stopMode = StopMode::CYC;
            else 
            {
                std::cerr << "Unknown stopMode: \"" << stopStr 
                          << "\". Defaulting to STEP.\n";
                stopMode = StopMode::STEP;
            }
        }

        if (((analysisMode == AnalysisMode::HB && numH != 0) || analysisMode == AnalysisMode::STEADY) &&
            stopMode != StopMode::STEP)
        {
            std::cout << "For analysisMode = HB or STEADY, stopMode must be STEP. "
                      << "Overriding user setting.\n";
            stopMode = StopMode::STEP;
        }

        if (config.count("writeMode")) {
            std::string writeStr = config["writeMode"];
            trim(writeStr);

            if (writeStr == "CYC")        writeMode = WriteMode::CYC;
            else if (writeStr == "STEP")   writeMode = WriteMode::STEP;
            else 
            {
                std::cerr << "Unknown writeMode: \"" << writeStr 
                          << "\". Defaulting to STEP.\n";
                writeMode = WriteMode::STEP;
            }
        }

        if (config.count("resWriteMode")) {
            std::string resWriteStr = config["resWriteMode"];
            trim(resWriteStr);

            if (resWriteStr == "CYC")        resWriteMode  = ResWriteMode::CYC;
            else if (resWriteStr == "STEP")   resWriteMode = ResWriteMode::STEP;
            else 
            {
                std::cerr << "Unknown resWriteMode: \"" << resWriteStr 
                          << "\". Defaulting to STEP.\n";
                resWriteMode = ResWriteMode::STEP;
            }
        }

        if (config.count("printMode")) {
            std::string printStr = config["printMode"];
            trim(printStr);

            if (printStr == "CYC")        printMode = PrintMode::CYC;
            else if (printStr == "STEP")   printMode = PrintMode::STEP;
            else 
            {
                std::cerr << "Unknown printMode: \"" << printStr 
                          << "\". Defaulting to STEP.\n";
                printMode = PrintMode::STEP;
            }
        }



        if (config.count("kB")) kB     = std::stod(config["kB"]);

        if (config.count("e")) e       = std::stod(config["e"]);

        if (config.count("eps0")) eps0 = std::stod(config["eps0"]);

        if (config.count("eV2J")) eV2J = std::stod(config["eV2J"]); 

        if (config.count("Hl")) Hl     = std::stod(config["Hl"]);

        Hl = Hl * eV2J;

        if (config.count("me")) me     = std::stod(config["me"]);

        if (config.count("Gam")) Gam   = std::stod(config["Gam"]);

        if (config.count("P")) P       = std::stod(config["P"]);

        if (config.count("N0")) N0     = std::stod(config["N0"]);

        if (config.count("De0")) De0   = std::stod(config["De0"]);

        if (config.count("Di0")) Di0   = std::stod(config["Di0"]);

        if (config.count("muE0")) muE0 = std::stod(config["muE0"]);

        if (config.count("muI0")) muI0 = std::stod(config["muI0"]);

        N = N0 * P;

        De = De0 / N;

        Di = Di0 / N;

        muE = muE0 / N;

        muI = muI0 / N;

        
        scalar theVelocity = std::sqrt(8.0 * kB * TeRef / (pi * me));

        tRef = LRef / theVelocity;

        fRef = 1.0 / tRef;

        EeRef  = 3.0/2.0 * kB * nRef * TeRef;

        ERef = phiRef / LRef;

        DeRef = LRef * theVelocity;

        DiRef = LRef * theVelocity;

        muERef = DeRef / phiRef;

        muIRef = DiRef / phiRef;

        eRef = EeRef /(phiRef * nRef);

        klRef = 1/(nRef * tRef);
        
        HlRef = EeRef / nRef;

        epsRef = EeRef * pow(LRef,2.0) / pow(phiRef,2.0);


        JeRef = nRef * theVelocity;

        JiRef = nRef * theVelocity;

        cSRef = klRef * pow(nRef,2.0);

        omegaRef = 2.0 * pi / tRef;

        if (config.count("stopTime")) stopTime           = std::stod(config["stopTime"]);

        if (config.count("stopStep")) stopStep           = std::stoi(config["stopStep"]);

        if (config.count("stopPeriod")) stopPeriod       = std::stoi(config["stopPeriod"]);

        if (config.count("writeCyc")) writeCyc           = std::stod(config["writeCyc"]);

        if (config.count("writeStep")) writeStep         = std::stoi(config["writeStep"]);

        if (config.count("resWriteCyc")) resWriteCyc     = std::stod(config["resWriteCyc"]);

        if (config.count("resWriteStep")) resWriteStep   = std::stoi(config["resWriteStep"]);

        if (config.count("printCyc")) printCyc           = std::stod(config["printCyc"]);

        if (config.count("printStep")) printStep         = std::stoi(config["printStep"]);

        if (config.count("nCellTs")) nCellTs             = std::stoi(config["nCellTs"]);

        if (config.count("cellIDTs")) {
            cellIDTs.clear();
            std::stringstream ss(config["cellIDTs"]);
            std::string val;
            while (std::getline(ss, val, ',')) {
                cellIDTs.push_back(std::stoi(val));
            }
        }

        if (config.count("nTI")) nTI = std::stoi(config["nTI"]);

        #define PRINT_VAR(v) \
            std::cout << std::left << std::setw(18) << #v << ": " << v << '\n';
        
        #define PRINT_LABEL(lab, v) \
        std::cout << std::left << std::setw(18) << lab << ": " << v << '\n';
        
        #define PRINT_ENUM(v) \
            std::cout << std::left << std::setw(18) << #v << ": " << toString(v) << '\n';

        #define PRINT_VEC(label, vec)                                            \
            do {                                                                 \
                std::cout << std::left << std::setw(18) << label << ": ";        \
                bool first__ = true;                                             \
                for (const auto &x__ : vec) {                                    \
                    if (!first__) std::cout << ' ';                              \
                    std::cout << x__;                                            \
                    first__ = false;                                             \
                }                                                                \
                std::cout << '\n';                                               \
            } while (0)

        std::cout << "================ Config Loaded ================" << std::endl;
        PRINT_VAR(initFile);
        PRINT_VAR(chemCSV);
        PRINT_VAR(meshFile);
        PRINT_VAR(numH);
        PRINT_VAR(numT);
        PRINT_VAR(PhiRF)
        PRINT_VAR(fRF)
        PRINT_LABEL("phase (radians)", phase)
        PRINT_VAR(period)
        PRINT_VAR(nRK)
        PRINT_VEC("alpha", alpha);
        PRINT_VAR(CFL)
        PRINT_VAR(eps)
        PRINT_VAR(nBJ)
        PRINT_VAR(tolNe)
        PRINT_VAR(tolNi)
        PRINT_VAR(tolEe)
        PRINT_VAR(Ni0)
        PRINT_VAR(Ne0)
        PRINT_VAR(Te0)
        PRINT_VAR(Phi0)
        PRINT_VAR(LRef)
        PRINT_VAR(TeRef)
        PRINT_VAR(tRef)
        PRINT_VAR(fRef)
        PRINT_VAR(EeRef)
        PRINT_VAR(nRef)
        PRINT_VAR(phiRef)
        PRINT_VAR(DeRef)
        PRINT_VAR(DiRef)
        PRINT_VAR(muERef)
        PRINT_VAR(muIRef)
        PRINT_VAR(eRef)
        PRINT_VAR(klRef)
        PRINT_VAR(HlRef)
        PRINT_VAR(epsRef)
        PRINT_VAR(cSRef)
        PRINT_VAR(JeRef)
        PRINT_VAR(JiRef)
        PRINT_VAR(ERef)
        PRINT_VAR(omegaRef)
        PRINT_VAR(isFirstOrder)
        PRINT_VAR(isJeConstraint)
        PRINT_ENUM(analysisMode)
        PRINT_ENUM(initOption)
        PRINT_ENUM(implicitScheme)
        PRINT_ENUM(timeStepType)
        PRINT_ENUM(stopMode)
        PRINT_ENUM(writeMode)
        PRINT_ENUM(resWriteMode)
        PRINT_ENUM(printMode)
        PRINT_VAR(pi)
        PRINT_VAR(kB)
        PRINT_VAR(e)
        PRINT_VAR(eps0)
        PRINT_VAR(eV2J)
        PRINT_VAR(Hl)
        PRINT_VAR(N)
        PRINT_VAR(De)
        PRINT_VAR(Di)
        PRINT_VAR(muE)
        PRINT_VAR(muI)
        PRINT_VAR(me)
        PRINT_VAR(Gam)
        PRINT_VAR(P)
        PRINT_VAR(N0)
        PRINT_VAR(De0)
        PRINT_VAR(Di0)
        PRINT_VAR(muE0)
        PRINT_VAR(muI0)
        PRINT_VAR(stopTime)
        PRINT_VAR(stopStep)
        PRINT_VAR(stopPeriod)
        PRINT_VAR(writeCyc)
        PRINT_VAR(writeStep)
        PRINT_VAR(resWriteCyc)
        PRINT_VAR(resWriteStep)
        PRINT_VAR(printCyc)
        PRINT_VAR(printStep)
        PRINT_VAR(nCellTs)
        PRINT_VEC("cellsIDTs", cellIDTs);
        PRINT_VAR(nTI)
        std::cout << "===============================================" << std::endl;
        std::cout << "Configuration loaded successfully!" << std::endl;        
    }



    void readMeshFile(const std::string &meshFile)
    {
        std::cout << "================ Mesh Loaded ================" << std::endl;
        std::ifstream infile(meshFile);
        if (!infile) {
            std::cerr << "Error opening mesh file: " << meshFile << std::endl;
            return;
        }

        std::string line;
        // Skip the first two comment lines
        std::getline(infile, line);
        std::getline(infile, line);

        // Read number of nodes
        std::getline(infile, line);
        numNodes = std::stoi(line);
        numCells = numNodes - 1;
        coords.resize(numNodes);

        // Read each coordinate line
        for (label i = 0; i < numNodes; ++i) {
            std::getline(infile, line);
            coords[i].x = std::stod(line);
        }
        infile.close();

        L = coords[numNodes - 1].x - coords[0].x;

        // Calculate the blocksize
        if (implicitScheme == ImplicitScheme::PCI1)
        {
            blockSizeS = 3;
            blockSizeM = blockSizeS * numCells;
            blockSizeL = blockSizeM * numT;
        }
        else if (implicitScheme == ImplicitScheme::FCI)
        {
            blockSizeS = 4;
            blockSizeM = blockSizeS * numCells;
            blockSizeL = blockSizeM * numT;
        } 
     
        std::cout << "Number of nodes: " << numNodes << std::endl;
        std::cout << "Number of cells: " << numCells << std::endl;
        std::cout << "Length of the domain: " << L << " m" << std::endl;
        std::cout << "Block sizes are: "
                  << blockSizeS << ", "
                  << blockSizeM << ", "
                  << blockSizeL << " respectively." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << "Mesh file loaded successfully!" << std::endl;

    }


    void scaleToDimensionless()
    {
        L      /= LRef;

        for (label i = 0; i < numNodes; ++i) {
            coords[i].x /= LRef;
        }

        e      /= eRef;

        eps0   /= epsRef;

        N      /= nRef;

        Hl     /= HlRef;

        De     /= DeRef;

        Di     /= DiRef;

        muE    /= muERef;

        muI    /= muIRef;

        PhiRF  /= phiRef;

        fRF    /= fRef;

        period /= tRef;

        Ni0    /= nRef;

        Ne0    /= nRef;

        Te0    /= TeRef;

        Phi0   /= phiRef;

        std::cout << "======= Dimensionless variables are as below =======" << std::endl;
        std::cout << "e: " << e << std::endl;
        std::cout << "eps0: " << eps0 << std::endl;
        std::cout << "N: " << N << std::endl;
        std::cout << "Hl: " << Hl << std::endl;
        std::cout << "De: " << De << std::endl;
        std::cout << "Di: " << Di << std::endl;
        std::cout << "muE: " << muE << std::endl;
        std::cout << "muI: " << muI << std::endl;
        std::cout << "PhiRF: " << PhiRF << std::endl;
        std::cout << "fRF: " << fRF << std::endl;
        std::cout << "period: " << period << std::endl;
        std::cout << "Ni0: " << Ni0 << std::endl;
        std::cout << "Ne0: " << Ne0 << std::endl;
        std::cout << "Te0: " << Te0 << std::endl;
        std::cout << "Phi0: " << Phi0 << std::endl;
        std::cout << "=====================================================" << std::endl;
        std::cout << "Nondimensionlization successfully done!" << std::endl;
    }

    void writeNodeLinePlot(const std::string& filename)
    {
        std::ofstream out(filename);
        if (!out.is_open())
        {
            std::cerr << "Error: Cannot open " << filename << " for writing node line plot.\n";
            return;
        }

        out << "TITLE = Node Distribution\n";
        out << "VARIABLES = X(m), Y\n";
        out << "ZONE T = Line, F=POINT\n";

        for (label i = 0; i < numNodes; ++i)
        {
            out << std::setprecision(8) << coords[i].x << "\t" << 0.0 << "\n";
        }

        out.close();
        std::cout << "Node line plot written to: " << filename << "\n";
    }

    


    std::string toString(AnalysisMode mode) {
        switch (mode) {
            case AnalysisMode::STEADY: return "STEADY";
            case AnalysisMode::DT:     return "DT";
            case AnalysisMode::HB:     return "HB";
            default:                   return "UNKNOWN";
        }
    }

    std::string toString(InitOption option) {
        switch (option) {
            case InitOption::SCRATCH:     return "SCRATCH";
            case InitOption::CONTINUE:    return "CONTINUE";
            case InitOption::RESTART:     return "RESTART";
            case InitOption::FROMSTEADY:  return "FROMSTEADY";
            default:                      return "UNKNOWN";
        }
    }

    std::string toString(ImplicitScheme scheme) {
        switch (scheme) {
            case ImplicitScheme::NO:   return "NO";
            case ImplicitScheme::PCI1: return "PCI1";
            case ImplicitScheme::FCI:  return "FCI";
            default:                   return "UNKNOWN";
        }
    }

    std::string toString(TimeStepType type) {
        switch (type) {
            case TimeStepType::LOCAL:  return "LOCAL";
            case TimeStepType::GLOBAL: return "GLOBAL";
            default:               return "UNKNOWN";
        }
    }

    std::string toString(StopMode mode) {
        switch (mode) {
            case StopMode::TIME:   return "TIME";
            case StopMode::STEP:   return "STEP";
            case StopMode::CYC:    return "CYC";
            default:               return "UNKNOWN";
        }
    }

    std::string toString(WriteMode mode) {
        switch (mode) {
            case WriteMode::CYC:   return "CYC";
            case WriteMode::STEP:  return "STEP";
            default:               return "UNKNOWN";
        }
    }

    std::string toString(ResWriteMode mode) {
        switch (mode) {
            case ResWriteMode::CYC:   return "CYC";
            case ResWriteMode::STEP:  return "STEP";
            default:               return "UNKNOWN";
        }
    }

    std::string toString(PrintMode mode) {
        switch (mode) {
            case PrintMode::CYC:   return "CYC";
            case PrintMode::STEP:  return "STEP";
            default:               return "UNKNOWN";
        }
    }


}