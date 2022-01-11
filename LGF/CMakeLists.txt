# CMakeList.txt : CMake project for LGF, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("LGF")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Add source to this project's executable.
add_executable (LGF "LGF.cpp" "LGF.h"   "Action/Action.h" "Action/Action.cpp" "Action/Gauge/Wilson.h" "Action/Gauge/Wilson.cpp" "Field/Field.h"  "Field/SU3_field.h" "Field/SU3_field.cpp" "Lattice/Lattice.h" "Lattice/Lattice.cpp" "Vertex.h" "Vertex.cpp" "mpiWrapper.h" "mpiWrapper.cpp"  "Math/SU3_mat.h" "Math/SU3_vec.h" "Math/SU3_mat.cpp" "Math/C_double.h" "Math/C_double.cpp" "Math/Random.h" "HMC/HMC.h" "HMC/HMC.cpp")
target_link_libraries(LGF ${MPI_LIBRARIES})
# TODO: Add tests and install targets if needed.
