// LGF.h : Include file for standard system include files,
// or project specific include files.
#pragma once

#include <iostream>
#include <climits>

// TODO: Reference additional headers your program requires here.
#include "Action/Action.h"
#include "Action/Gauge/Wilson.h"
#include "Field/Field.h"
#include "Field/SU3_field.h"
#include "Lattice/Lattice.h"
#include "mpiWrapper.h"
#include "mpi.h"
#include "Math/Random.h"
#include "HMC/HMC.h"
#include "HMC/LHMC.h"
#include "LMC/LMC.h"
#include "Math/SU3_gen.h"
#include "Gradient_Flow/GradientFlow.h"
#include "Observables/TopologicalCharge.h"
#include "Observables/EnergyDensity.h"
#include "Observables/QdensityDistribution.h"
#include "SU3_Heatbath/SU3_Heatbath.h"

void GenerateHeatBathGaugeConfigurations(SimulationParameters& params, HBParameters& HBparams);
void FlowSavedGaugeConfigurations(SimulationParameters& params, FlowParameters& flowParams);
void ResumeFlowedConfiguration(SimulationParameters& params, FlowParameters& flowParams, double flowTime_pickup);

//void MeasureFlowedGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, double flowTime);
//void ResumeHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationResumeFrom, int ConfigurationStop, int ConfigTimeSeparation);
//#define mpi_debug_breakpoint int temp;if (mpiWrapper::id() == 0) {std::cout << "input number to continue application: ";std::cin >> temp;} MPI_Barrier(mpiWrapper::comm());

void ParseCLargs(int argc, char** argv, int& NrDims, int& extdofs, double& beta, std::vector<int>& shape, std::vector<int>& cuts, int& ConfigStart, int& ConfigStop, int& flowSteps, double& epsilon, int& measureInterval, int& ThermSteps, int& OR_per_HB, int& configSep, std::string& dataFolder);