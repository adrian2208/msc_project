// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include "filesystem.h"


int main(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	//Physical Parameters
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = {16, 16, 16, 16};
	double beta = 6.0;
	int cuts[] = { 1,1,1,0 }; // shape_x must be divisible by (cut_x + 1)
	SimulationParameters params(NrDims, extdofs, shape, beta, cuts);
	int ConfigurationStart = 0;
	int ConfigurationStop = 300;

	//FlowSavedGaugeConfigurations(params, ConfigurationStart, ConfigurationStop,350,5);//300

	int ThermalizationSteps = 2000;
	int ConfigTimeSeparation = 200;
	//GenerateHeatBathGaugeConfigurations(params, ConfigurationStart, ConfigurationStop, ThermalizationSteps, ConfigTimeSeparation);
	int ConfigurationResumeFrom = 267;
	ResumeHeatBathGaugeConfigurations(params, ConfigurationResumeFrom, ConfigurationStop,ConfigTimeSeparation);


	//ResumeFlowedConfiguration(params, ConfigurationStart, ConfigurationStop, 200,16.0);

	double measureAtFlowTime = 8.02;
	//MeasureFlowedGaugeConfigurations(params, ConfigurationStart, ConfigurationStop, measureAtFlowTime);

	//MEASURE Q_DENSITYDISTRIBUTION
	//Lattice lattice(NrDims, shape);
	//SU3_field U(lattice, extdofs);
	//std::string ensembleNum = std::to_string(170);
	//double flowTime = 20.0;
	//U.loadSU3FromFile(beta, "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime));
	//QdensityDistribution Qdensity(U, beta, "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime));
	//Qdensity.calculate(flowTime);

	mpiWrapper::end_parallelSession();

	return 0;
}

void MeasureFlowedGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, double flowTime) {
	Lattice lattice(params, true);
	Wilson action(params.getbeta());
	//run through the saved gauge configurations
	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime));

		//instantiate and include the Topological Charge
		TopologicalCharge topCharge(U);
		//topCharge.calculate(0.0);
		//instantiate and include the Energy density
		//EnergyDensity Edensity(U);

		//save the observables to files
		topCharge.saveTopologicalChargeToFile(params.getbeta(),"GF" , ensembleNum + "_Flowtime" + std::to_string(flowTime));
		//Edensity.saveEnergyDensityToFile(beta, ensembleNum + "flowed");
	}
}

void ResumeHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationResumeFrom, int ConfigurationStop,int ConfigTimeSeparation) {
	Lattice lattice(params,false);
	SU3_field U(lattice, params.getextdofs());
	SU3_Heatbath heatbath(U, params.getbeta(), 4);
	U.loadSU3FromFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(ConfigurationResumeFrom));
	for (int i = ConfigurationResumeFrom+1; i <= ConfigurationStop; i++) {
		heatbath.update(ConfigTimeSeparation);
		U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(i));
	}
}
void GenerateHeatBathGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, int ThermalizationSteps, int ConfigTimeSeparation) {
	Lattice lattice(params,false);
	SU3_field U(lattice, params.getextdofs());
	U.InitializeHotStart();
	SU3_Heatbath heatbath(U, params.getbeta(), 4);
	heatbath.update(ThermalizationSteps);
	U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(ConfigurationStart));
	for (int i = ConfigurationStart+1; i <= ConfigurationStop; i++) {
		heatbath.update(ConfigTimeSeparation);
		U.saveSU3ToFile(params.getbeta(), heatbath.getupdateMethod(), std::to_string(i));
	}
}
void FlowSavedGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, int flowSteps, int measure_every_nth_step) {
	//flow step size
	double epsilon = 0.01;
	//instantiate the lattice and action
	Lattice lattice(params,true);
	Wilson action(params.getbeta());
	//run through the saved gauge configurations
	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "Heatbath_4_ORperHB", ensembleNum);
		//instantiate the Gradient flow object
		GradientFlow flowing(action, U, epsilon, measure_every_nth_step);

		
		//instantiate and include the Topological Charge
		TopologicalCharge topCharge(U);
		flowing.Include_TopCharge(topCharge);
		//instantiate and include the Energy density
		EnergyDensity Edensity(U);
		flowing.Include_EnergyDensity(Edensity);


		double flowTime;
		int SaveEveryNconfigs = 500;
		//flow the configuration
		for (int i = 0; i < flowSteps; i++) {
			flowing.flow();
			flowTime = flowing.GetFlowTime();
		}
		////save the observables to files
		topCharge.saveTopologicalChargeToFile(params.getbeta(),flowing.getupdateMethod() ,"_" + ensembleNum);
		Edensity.saveEnergyDensityToFile(params.getbeta(), flowing.getupdateMethod(), "_" + ensembleNum);
		U.saveSU3ToFile(params.getbeta(), flowing.getupdateMethod(), ensembleNum + "_Flowtime" + std::to_string(flowTime));
	}
}
void ResumeFlowedConfiguration(SimulationParameters& params, int ConfigurationStart,int ConfigurationStop, int extra_updates,double flowTime) {
	//flow step size
	double epsilon = 0.02;
	//instantiate the lattice and action
	Lattice lattice(params, true);
	Wilson action(params.getbeta());
	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, params.getextdofs());
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(params.getbeta(), "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime));
		//instantiate the Gradient flow object
		GradientFlow flowing(action, U, epsilon);
		//instantiate and include the Topological Charge
		//TopologicalCharge topCharge(U);
		//flowing.Include_TopCharge(topCharge);

		//flow the configuration
		for (int i = 0; i < extra_updates; i++) {
			flowing.flow();
		}
		//save the observables to files
		//topCharge.saveTopologicalChargeToFile(beta, flowing.getupdateMethod(), "Heatbath_4_ORperHB_Flowed" + ensembleNum);
		//Edensity.saveEnergyDensityToFile(beta, ensembleNum + "flowed");
		double NewflowTime = flowing.GetFlowTime();
		U.saveSU3ToFile(params.getbeta(), "GF", ensembleNum + "_Flowtime" + std::to_string(flowTime + NewflowTime));
	}
}

void GenerateLHMCGaugeConfigurations(SimulationParameters& params, int ConfigurationStart, int ConfigurationStop, int NrUpdates) {
	
	double LHMCepsilon = 0.33;
	int NrLeaps = 3;

	Lattice lattice(params, false);
	Wilson action(params.getbeta());

	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		std::string ensembleNum = "ensemble_"+ std::to_string(i);
		SU3_field U(lattice, params.getextdofs());
		U.InitializeHotStart();
		LHMC updater(U, action, LHMCepsilon,NrLeaps);

		TopologicalCharge topCharge(U);
		EnergyDensity Edensity(U);
		for (int i = 0; i < NrUpdates; i++) {
			updater.sweep();
			U.transfer_FieldValues();
			topCharge.calculate(0);
			Edensity.calculate(0);
		}
		std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";
		MPI_Barrier(mpiWrapper::comm());
		U.saveSU3ToFile(params.getbeta(), ensembleNum);
		topCharge.saveTopologicalChargeToFile(params.getbeta(), ensembleNum +"LHMC");
		Edensity.saveEnergyDensityToFile(params.getbeta(), ensembleNum+"LHMC");
	}
}


