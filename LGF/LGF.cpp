#include "LGF.h"
#include <cstring>
#ifdef LEGACY_CXX
#include <experimental/filesystem>
namespace n_fs = ::std::experimental::filesystem;
#else
#include <filesystem>
namespace n_fs = ::std::filesystem;
#endif


//TEMP
#include <chrono>
#include <thread>

using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
using std::chrono::system_clock;

/*
 RUN EXAMPLE:
	*edit main scope to include either
		FlowSavedGaugeConfigurations(params, flowParams);
						or
		GenerateHeatBathGaugeConfigurations(params, HBparams);
		
		or some other, less often used function
	*compile

	*WITH COMMAND LINE ARGS:
	mpiexec.exe -n Nprocs LGF.exe NrDims extdofs beta shape cuts configStart configStop flowSteps/ThermSteps epsilon/OR_perHB MeasureInterval/configSep dataFolder
	
	*WITHOUT COMMAND LINE ARGS:
	mpiexec.exe -n Nprocs LGF.exe
		--- edit parameters before compilation
 
 */
int main(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	//Physical Parameters
	int NrDims; int extdofs; double beta; std::vector<int> shape;
	//Partition Parameter
	std::vector<int> cuts; // shape_x must be divisible by (cut_x + 1)
	//Flow Parameters
	int ConfigStart; int ConfigStop; int flowSteps; double epsilon; int measureInterval;
	//Heatbath Parameters
	int ThermSteps; int OR_per_HB; int configSep;
	//output folder for observables and configurations
	std::string dataFolder;

	if (argc > 1) {
		ParseCLargs(argc, argv, NrDims, extdofs, beta, shape, cuts, ConfigStart, ConfigStop, flowSteps, epsilon, measureInterval, ThermSteps, OR_per_HB, configSep, dataFolder);
	}
	else {
		////Physical Parameters
		//NrDims = 4;
		//extdofs = 4;
		//beta = 6.26;
		//shape = { 24,24,24,24 };
		////Partition Parameter
		//cuts = { 2,1,1,0 }; // shape_x must be divisible by (cut_x + 1)
		////File Parameters
		//ConfigStart =0;
		//ConfigStop = 0;
		////Flow Parameters
		//flowSteps = 10;
		//epsilon = 0.01;
		//measureInterval = 10;
		////Heatbath Parameters
		//ThermSteps = 1000;
		//OR_per_HB = 4;
		//configSep = 1;
		////output folder for observables and configurations
		//dataFolder = "C:/Users/adria/Documents/msc_project/data/";
				//Physical Parameters
		NrDims = 4;
		extdofs = 4;
		beta = 6.13;
		shape = { 12,12,12,12 };
		//Partition Parameter
		cuts = { 3,0,0,0 }; // shape_x must be divisible by (cut_x + 1)
		//File Parameters
		ConfigStart = 0;
		ConfigStop = 100;
		//Flow Parameters
		flowSteps = 40;
		epsilon = 0.33;
		measureInterval = 1;
		//Heatbath Parameters
		ThermSteps = 1000;
		OR_per_HB = 4;
		configSep = 1;
		//output folder for observables and configurations
		dataFolder = "C:/Users/adria/Documents/msc_project/data/";
	}
	SimulationParameters params(NrDims, extdofs, shape, beta, cuts,"torus");//_OpenT
	FlowParameters flowParams(ConfigStart, ConfigStop, flowSteps, epsilon, measureInterval, dataFolder);
	HBParameters HBparams(ConfigStart, ConfigStop, ThermSteps, OR_per_HB, configSep, dataFolder);
	LHMCParameters LHMCparams(ConfigStart, ConfigStop, ThermSteps, epsilon, configSep, dataFolder);
	//FlowSavedGaugeConfigurations(params, flowParams);

	double flowTime_pickup = 20;
	//ResumeFlowedConfiguration(params, flowParams, flowTime_pickup);
	//MeasureFlowedGaugeConfigurations(params, flowParams);
	//MeasureQdensityDistribution(params, flowParams);
	//GenerateHeatBathGaugeConfigurations(params, HBparams);
	GenerateLHMCGaugeConfigurations(params, LHMCparams);

	mpiWrapper::end_parallelSession();

	return 0;
}

void MeasureQdensityDistribution(SimulationParameters& params, FlowParameters& flowParams) {
	Lattice lattice(params,true);
	double flowTime_pickup = flowParams.getepsilon() * flowParams.getflowSteps();
	for (int i = flowParams.getConfigurationStart(); i <= flowParams.getConfigurationStop(); i++) {
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "GF", flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(flowTime_pickup));
		QdensityDistribution Qdensity(U, params.getbeta(), "GF", flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(flowTime_pickup));
		Qdensity.calculate(flowTime_pickup);
	}
}

//void ResumeHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationResumeFrom, int ConfigurationStop,int ConfigTimeSeparation) {
//	Lattice lattice(params,false);
//	SU3_field U(lattice, params.getextdofs());
//	SU3_Heatbath heatbath(U, params.getbeta(), 4);
//	U.loadSU3FromFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(ConfigurationResumeFrom));
//	for (int i = ConfigurationResumeFrom+1; i <= ConfigurationStop; i++) {
//		heatbath.update(ConfigTimeSeparation);
//		U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(i));
//	}
//}
void GenerateHeatBathGaugeConfigurations(SimulationParameters& params, HBParameters& HBparams) {
	Lattice lattice(params,false);
	SU3_field U(lattice, params.getextdofs());
	//Check whether this is the start of a new ensemble, which needs thermalization
	bool NewEnsemble = HBparams.getThermSteps() != 0 && HBparams.getConfigurationStart() == 0;
	//If so, initialize from a random configuration
	if (NewEnsemble) {
		U.InitializeHotStart();
	}
	//otherwise, load the previously generated configuration to continue 
	else {
		U.loadSU3FromFile(params.getbeta(), "Heatbath_" + std::to_string(HBparams.getOR_per_HB()) + "_ORperHB", HBparams.getdataFolder(), std::to_string(HBparams.getConfigurationStart()));
	}
	//initialize the update program
	SU3_Heatbath heatbath(U, params.getbeta(), HBparams.getOR_per_HB());
	//update the field either for thermalization or configuration separation depending on what the user is doing
	if (NewEnsemble) {
		heatbath.update(HBparams.getThermSteps());
		U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), HBparams.getdataFolder(), std::to_string(HBparams.getConfigurationStart()));
	}
	for (int i = HBparams.getConfigurationStart() + 1; i <= HBparams.getConfigurationStop(); i++) {
		heatbath.update(HBparams.getconfigSep());
		U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), HBparams.getdataFolder(), std::to_string(i));
	}
}


void FlowSavedGaugeConfigurations(SimulationParameters& params, FlowParameters& flowParams) {
	std::string UpdateMethod = "LHMC_3_of_0_330000";
	//instantiate the lattice and action
	Lattice lattice(params,true);
	Wilson action(params.getbeta());
	//run through the saved gauge configurations
	for (int i = flowParams.getConfigurationStart(); i <= flowParams.getConfigurationStop(); i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), UpdateMethod, flowParams.getdataFolder(), ensembleNum);
		//instantiate the Gradient flow object
		GradientFlow flowing(action, U, flowParams.getepsilon(), flowParams.getMeasuringInterval());
		
		//instantiate and include the Topological Charge
		TopologicalCharge topCharge(U);
		flowing.Include_TopCharge(topCharge);
		topCharge.calculate(0.0);
		//instantiate and include the Energy density
		EnergyDensity Edensity(U);
		flowing.Include_EnergyDensity(Edensity);
		Edensity.calculate(0.0);
		//flow the configuration
		for (int i = 0; i < flowParams.getflowSteps(); i++) {
			flowing.flow();
		}
		double flowTime = flowing.GetFlowTime();
		//save the observables to files
		topCharge.saveTopologicalChargeToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(),"_" + ensembleNum);
		Edensity.saveEnergyDensityToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(), "_" + ensembleNum);
		//Save the last configuration to file
		U.saveSU3ToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(flowTime));
	}
}
void ResumeFlowedConfiguration(SimulationParameters& params, FlowParameters& flowParams,double flowTime_pickup) {
	//instantiate the lattice and action
	Lattice lattice(params, true);
	Wilson action(params.getbeta());
	for (int i = flowParams.getConfigurationStart(); i <= flowParams.getConfigurationStop(); i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "GF", flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(flowTime_pickup));
		//instantiate the Gradient flow object
		GradientFlow flowing(action, U, flowParams.getepsilon(), flowParams.getMeasuringInterval());
		flowing.SetFlowTime(flowTime_pickup);
		//instantiate and include the Topological Charge
		TopologicalCharge topCharge(U);
		flowing.Include_TopCharge(topCharge);
		//instantiate and include the Energy density
		EnergyDensity Edensity(U);
		flowing.Include_EnergyDensity(Edensity);

		//flow the configuration
		for (int i = 0; i < flowParams.getflowSteps(); i++) {
			flowing.flow();
		}
		//save the observables to files
		topCharge.saveTopologicalChargeToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(), "_" + ensembleNum);
		Edensity.saveEnergyDensityToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(), "_" + ensembleNum);
		double NewflowTime = flowing.GetFlowTime();
		U.saveSU3ToFile(params.getbeta(), flowing.getupdateMethod(), flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(NewflowTime));
	}
}

void MeasureFlowedGaugeConfigurations(SimulationParameters& params, FlowParameters& flowParams) {
	Lattice lattice(params, true);
	Wilson action(params.getbeta());
	double flowTime_pickup = flowParams.getepsilon() * flowParams.getflowSteps();
	//run through the saved gauge configurations
	for (int i = flowParams.getConfigurationStart(); i <= flowParams.getConfigurationStop(); i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "GF", flowParams.getdataFolder(), ensembleNum + "_Flowtime" + std::to_string(flowTime_pickup));

		//instantiate and include the Topological Charge
		UnImprovedTopologicalCharge unimprovedTopCharge(U);
		unimprovedTopCharge.calculate(flowTime_pickup);
		//TopologicalCharge topCharge(U);
		//topCharge.calculate(0.0);
		//instantiate and include the Energy density
		//EnergyDensity Edensity(U);
		ImprovedEnergyDensity ImprovedEdensity(U);
		ImprovedEdensity.calculate(flowTime_pickup);

		//save the observables to files
		//unimprovedTopCharge.saveTopologicalChargeToFile(params.getbeta(), "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime));
		unimprovedTopCharge.saveTopologicalChargeToFile(params.getbeta(), "GF", flowParams.getdataFolder(), "_" + ensembleNum);
		ImprovedEdensity.saveEnergyDensityToFile(params.getbeta(), "GF", flowParams.getdataFolder(), "_" + ensembleNum);
		//Edensity.saveEnergyDensityToFile(beta, ensembleNum + "flowed");
	}
}

void GenerateLHMCGaugeConfigurations(SimulationParameters& params, LHMCParameters LHMCparams) {
	double LHMCepsilon = LHMCparams.getEpsilon();
	int NrLeaps = LHMCparams.getLeapfrogSteps();
	//sleep_for(30s);
	mpi_debug_breakpoint
	Lattice lattice(params, false);
	Wilson action(params.getbeta());
	SU3_field U(lattice, params.getextdofs());
	//initialize the update program
	LHMC updater(U, action, LHMCepsilon, NrLeaps);


	std::cout << "ProcID: " << mpiWrapper::id() << "		proc Volume: " << lattice.getthisProc_Volume() << "		resp Volume: " << lattice.getthisProc_RespVolume() << std::endl;
	lattice.Print_Partitioning();


	//Check whether this is the start of a new ensemble, which needs thermalization
	bool NewEnsemble = LHMCparams.getThermSteps() != 0 && LHMCparams.getConfigurationStart() == 0;

	//If so, initialize from a random configuration
	if (NewEnsemble) {
		U.InitializeHotStart();
	}
	//otherwise, load the previously generated configuration to continue 
	else {
		U.loadSU3FromFile(params.getbeta(), updater.getupdateMethod(), LHMCparams.getdataFolder(), std::to_string(LHMCparams.getConfigurationStart()));
	}

	//update the field either for thermalization or configuration separation depending on what the user is doing
	if (NewEnsemble) {
		for (int i = 0; i < LHMCparams.getThermSteps(); i++) {
			if (mpiWrapper::id() == 0) {
				std::cout << "step: " << i << std::endl;
			}
			
			updater.sweep();
		}
		U.saveSU3ToFile(params.getbeta(), updater.getupdateMethod(), LHMCparams.getdataFolder(), std::to_string(LHMCparams.getConfigurationStart()));
	}
	for (int i = LHMCparams.getConfigurationStart() + 1; i <= LHMCparams.getConfigurationStop(); i++) {
		for (int j = 0; j < LHMCparams.getconfigSep(); j++) {
			if (mpiWrapper::id() == 0) {
				std::cout << "step: " << i << std::endl;
			}
			updater.sweep();
		}
		U.saveSU3ToFile(params.getbeta(), updater.getupdateMethod(), LHMCparams.getdataFolder(), std::to_string(i));
	}
}


void ParseCLargs(int argc, char** argv, int& NrDims, int& extdofs, double& beta, std::vector<int>& shape, std::vector<int>& cuts, int& ConfigStart, int& ConfigStop, int& flowSteps, double& epsilon, int& measureInterval, int& ThermSteps, int& OR_per_HB, int& configSep, std::string& dataFolder) {
	double Flow_v_LHMC_eps_Cutoff = 0.1;
	NrDims = std::stoi(argv[1]);
	extdofs = std::stoi(argv[2]);
	beta = std::stod(argv[3]);
	int i = 4;
	for (; i < 4 + NrDims; i++) {
		shape.push_back(std::stoi(argv[i]));
		cuts.push_back(std::stoi(argv[i + NrDims]));
	}
	i += 4;
	ConfigStart = std::stoi(argv[i]);
	i += 1;
	ConfigStop = std::stoi(argv[i]);
	i += 1;
	flowSteps = std::stoi(argv[i]);
	ThermSteps = std::stoi(argv[i]);
	i += 1;
	std::size_t temp;
	std::stoi(argv[i], &temp);
	if (temp == std::strlen(argv[i])) {
		OR_per_HB = std::stoi(argv[i]);
		epsilon = INT_MAX;
	}
	else {
		OR_per_HB = INT_MAX;
		epsilon = std::stod(argv[i]);
	}
	i += 1;
	measureInterval = std::stoi(argv[i]);
	configSep = std::stoi(argv[i]);
	i += 1;
	dataFolder = argv[i];

	if (mpiWrapper::id() == 0) {
		std::cout << "Command Line Args parsed as:\n" << "NrDims = " << NrDims << "\nextdofs = " << extdofs << "\nbeta = " << beta << "\nshape = ";
		for (auto& item : shape) {
			std::cout << std::to_string(item) << " ";
		}
		std::cout << "\ncuts = ";
		for (auto& item : cuts) {
			std::cout << std::to_string(item) << " ";
		}
		std::string ThermStatus;
		if (ThermSteps != 0) {
			ThermStatus = std::to_string(ThermSteps);
		}
		else {
			ThermStatus = "None";
		}
		std::cout << "\nconfigs: " + std::to_string(ConfigStart) + " --> " + std::to_string(ConfigStop) + "\n";
		if (OR_per_HB != INT_MAX) {
			std::cout << "Thermalization: " + ThermStatus + "	OR-steps per HB_step: " << std::to_string(OR_per_HB) << "	Config Separation: " << std::to_string(configSep) << "\n";
		}
		else {
			if (epsilon < Flow_v_LHMC_eps_Cutoff) {
				std::cout << "Epsilon: " << std::to_string(epsilon) << "	Measuring interval: " << std::to_string(measureInterval) << "	flowSteps: " << std::to_string(flowSteps) << "\n";
			}
			else {
				std::cout << "Thermalization: " + ThermStatus + "LHMC_Epsilon: " << std::to_string(epsilon) << "	Config Separation: " << std::to_string(configSep) << "\n";
			}
			
		}
		std::cout << "Data Directory: " + dataFolder << std::endl;
	}

}