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
#include "Vertex.h"
#include "mpiWrapper.h"
#include "mpi.h"
#include "Math/Random.h"
#include "HMC/HMC.h"
#include "LMC/LMC.h"
#include "Math/SU3_gen.h"



void testSaveLoad(int argc, char** argv);
void testHMC(int argc, char** argv);
void testLMC(int argc, char** argv);


#define mpi_debug_breakpoint int temp;if (mpiWrapper::id() == 0) {std::cout << "input number to continue application: ";std::cin >> temp;} MPI_Barrier(mpiWrapper::comm());