# Lattice Gradient Flow (LGF)
## A C++ tool for generating pure SU(3) gauge configurations, applying the Wilson flow to those configuration, and analyzing the energy density and topological charge at all flow-times

This program is accompanied by a thesis analyzing several aspects of the Wilson flow, the topological charge- and susceptibility, and the Energy density of pure SU(3) gauge configurations.
The thesis can be found in the 'doc'- directory.
## Building LGF

The LGF program contains several primary functions contained in LGF.cpp. These contain the primary functionalities of the program.
They are:
```c++
GenerateHeatBathGaugeConfigurations(params, HBparams);
FlowSavedGaugeConfigurations(params, flowParams);
ResumeFlowedConfiguration(params, flowParams, flowTime_pickup); //flowTime_pickup is not supported by command line arguments and must be compiled
MeasureQdensityDistribution(params, flowParams);
```

In order to run either of these four, all functions except the relevant function must be commented out. 
All of these are contained in the main scope between the declerations 
```c++
HBParameters HBparams(ConfigStart, ConfigStop, ThermSteps, OR_per_HB, configSep, dataFolder);
```
and 
```c++
mpiWrapper::end_parallelSession();
```
You may build the program in the following way:
```Shell
cd LGF
mkdir f_build ### Comment out all but FlowSavedGaugeConfigurations
mkdir g_build ### Comment out all but GenerateHeatBathGaugeConfigurations
mkdir r_build ### Comment out all but ResumeFlowedConfiguration
mkdir m_build ### Comment out all but MeasureQdensityDistribution

cd x_build
cmake --CMAKE_BUILD_TYPE=Release
cmake --build . ### You may have to edit set(CMAKE_CXX_STANDARD 20) in CMakeLists depending on your C++ standard
```
## Running LGF
LGF takes the same amount of command line arguments regardless of build. 
The exception is r_build, which doesn't currently support a command line argument for flowTime_pickup.
The command line arguments are entered in the following order
```Shell
mpiexec.exe -n Nprocs LGF.exe  NrDims  extdofs  beta  shape  cuts  configStart  configStop  flowSteps/ThermSteps  epsilon/OR_perHB  MeasureInterval/configSep  dataFolder
```
There is an interplay between Nprocs, shape and cuts which needs to be fulfilled in order for the program to execute succesfully.
In order to determine the appropriate values for 'cuts' and 'Nprocs', one can use the partitioning_check utility provided in the 'Simulation- Tools' directory.
