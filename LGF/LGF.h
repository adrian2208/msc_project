﻿// LGF.h : Include file for standard system include files,
// or project specific include files.
#pragma once

#include <iostream>

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

void GenerateHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, int ThermalizationSteps, int ConfigTimeSeparation);
void ResumeHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationResumeFrom, int ConfigurationStop, int ConfigTimeSeparation);
void FlowSavedGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, int flowSteps, int measure_every_nth_step);
void ResumeFlowedConfiguration(SimulationParameters& params, int ConfigurationStart,int ConfigurationStop, int extra_updates,double flowTime);
void MeasureFlowedGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, double flowTime);

//#define mpi_debug_breakpoint int temp;if (mpiWrapper::id() == 0) {std::cout << "input number to continue application: ";std::cin >> temp;} MPI_Barrier(mpiWrapper::comm());

