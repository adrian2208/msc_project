# CMakeList.txt : CMake project for LGF, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("LGF")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# EXPERIMENTING -- FIRST TO GO IN CASE OF FAILURE -----------<>
#set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -O3 -funroll-loops -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
# EXPERIMENTING -- FIRST TO GO IN CASE OF FAILURE -------------^

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Add source to this project's executable.
add_executable (LGF "LGF.cpp" "LGF.h"   "Action/Action.h" "Action/Action.cpp" "Action/Gauge/Wilson.h" "Action/Gauge/Wilson.cpp" "Field/Field.h"  "Field/SU3_field.h" "Field/SU3_field.cpp" "Lattice/Lattice.h" "Lattice/Lattice.cpp"   "mpiWrapper.h" "mpiWrapper.cpp"  "Math/SU3_mat.h" "Math/SU3_vec.h" "Math/SU3_mat.cpp" "Math/C_double.h" "Math/C_double.cpp" "Math/Random.h" "HMC/HMC.h" "HMC/HMC.cpp" "Field/RNG_field.h" "Field/RNG_field.cpp" "Math/SU3_gen.h" "Math/SU3_gen.cpp" "LMC/LMC.h" "LMC/LMC.cpp" "HMC/LHMC.cpp" "HMC/LHMC.h" "Field/Double_field.h" "Field/Double_field.cpp" "Gradient_Flow/GradientFlow.h" "Gradient_Flow/GradientFlow.cpp" "Observables/TopologicalCharge.cpp" "Observables/TopologicalCharge.h"  "SU3_Heatbath/SU3_Heatbath.h" "SU3_Heatbath/SU3_Heatbath.cpp" "Observables/EnergyDensity.h" "Observables/EnergyDensity.cpp" "Observables/ImprovedTopologicalCharge.h" "Observables/ImprovedTopologicalCharge.cpp" "Gradient_Flow/FastGF.cpp" "Observables/QdensityDistribution.h" "Observables/QdensityDistribution.cpp")
target_link_libraries(LGF ${MPI_LIBRARIES})

