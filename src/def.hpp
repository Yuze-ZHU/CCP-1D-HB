// Definition of the Globals class - GLOBAL variables and functions.

#pragma once
#include <cmath>
#include <cfloat>
#include <cstddef>
#include <string>

using label = long long int;
using scalar = double;

// General structures of geometry and of flow variables =======================


// Options for initializing the solution:
// SCRATCH     : Start from Ni0, Ne0, Te0, Phi0 (cold start)
// CONTINUE : Continue from a saved file, preserving the existing residual history
// RESTART  : Restart from a saved file, but reset and start a new residual history
// FROMSTEADY  : Initialize from a steady-state solution for unsteady simulation
enum class InitOption { SCRATCH, CONTINUE, RESTART, FROMSTEADY };

// Types of analysis modes:
// STEADY : Steady-state simulation (time-independent)
// DT     : Time-domain simulation (first-order backward Euler)
// HB     : Harmonic balance method (time-frequency-domain for periodic problems)
enum class AnalysisMode { STEADY, DT, HB };

// Implicit scheme for pseudo-time stepping:
// NO  	: No implicit scheme, all equations are solved explicitly
// FCI  ：Implicitly couples all governing equation including the Poisson equation
// PCI1 : Implicitly couples all governing equations except the Poisson equation
// PCI2 : Implicitly couples all gocerning equations except the Electric Energy equation
enum class ImplicitScheme { NO, FCI, PCI1, PCI2};

// Conditions to stop the simulation:
// TIME   : Stop when a specified physical time is reached
// STEP   : Stop after a certain number of time steps
// CYC    : Stop after a number of physical periods
enum class StopMode{ TIME, STEP, CYC };

// Type of time step used in the simulation:
// LOCAL  : Local time step (each cell may have a different step size)
// GLOBAL : Global time step (uniform for all cells)
enum class TimeStepType { LOCAL, GLOBAL };

// Type of output-write mode used in the simulation:
//
// CYC       : Write files every fixed number of running cycle
//             (controlled by writeCyc)
// STEP      : Write files every fixed number of running step
//             (controlled by writeStep)
enum class WriteMode { CYC, STEP };

// Type of residual output-write mode used in the simulation:
//
// CYC       : Write files every fixed number of running cycle
//             (controlled by writeCyc)
// STEP      : Write files every fixed number of running step
//             (controlled by writeStep)
enum class ResWriteMode { CYC, STEP };

// Type of print mode used in the simulation:
//
// CYC       : Print information on screen every fixed number of running cycle
//             (controlled by printCyc)
// STEP      : Print information on screen every fixed number of running step
//             (controlled by printStep)
enum class PrintMode { CYC, STEP };

// Coordinates of a grid node(1D version)
typedef struct T_NODE
{
	scalar x; // x-coordinate
} NODE;

// Components of a face-normal vector in 1D.
typedef struct T_VECTOR
{
	double nx;
} VECTOR;
