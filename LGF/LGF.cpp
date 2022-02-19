// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>


int main(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	//Physical Parameters
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = {24,24,24,24};
	double beta = 6.0;
	int ConfigurationStart = 81;
	int ConfigurationStop = 100;
	//std::cout << getLatticeConstant(beta)<< std::endl;
	//GenerateLHMCGaugeConfigurations(NrDims, extdofs, shape, beta, ConfigurationStart, ConfigurationStop, 200);
	FlowSavedGaugeConfigurations(NrDims, extdofs, shape, beta, ConfigurationStart, ConfigurationStop,800);

	//int ThermalizationSteps = 0;
	int ConfigTimeSeparation = 200;
	//GenerateHeatBathGaugeConfigurations(NrDims, extdofs, shape, beta, ConfigurationStart, ConfigurationStop, ThermalizationSteps, ConfigTimeSeparation);
	int ConfigurationResumeFrom = 80;
	//ResumeHeatBathGaugeConfigurations(NrDims, extdofs, shape, beta, ConfigurationResumeFrom, ConfigurationStop,ConfigTimeSeparation);

	mpiWrapper::end_parallelSession();

	return 0;
}
void ResumeHeatBathGaugeConfigurations(int NrDims, int extdofs, int shape[], double beta, int ConfigurationResumeFrom, int ConfigurationStop,int ConfigTimeSeparation) {
	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	SU3_Heatbath heatbath(U, beta, 4);
	U.loadSU3FromFile(beta, heatbath.getupdateMethod(), std::to_string(ConfigurationResumeFrom));
	for (int i = ConfigurationResumeFrom+1; i <= ConfigurationStop; i++) {
		heatbath.update(ConfigTimeSeparation);
		U.saveSU3ToFile(beta, heatbath.getupdateMethod(), std::to_string(i));
	}
}
void GenerateHeatBathGaugeConfigurations(int NrDims, int extdofs, int shape[], double beta, int ConfigurationStart, int ConfigurationStop, int ThermalizationSteps, int ConfigTimeSeparation) {
	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	SU3_Heatbath heatbath(U, beta, 4);
	heatbath.update(ThermalizationSteps);
	U.saveSU3ToFile(beta, heatbath.getupdateMethod(), std::to_string(ConfigurationStart));
	for (int i = ConfigurationStart+1; i <= ConfigurationStop; i++) {
		heatbath.update(ConfigTimeSeparation);
		U.saveSU3ToFile(beta, heatbath.getupdateMethod(), std::to_string(i));
	}
}
void FlowSavedGaugeConfigurations(int NrDims, int extdofs, int shape[], double beta, int ConfigurationStart, int ConfigurationStop, int flowSteps) {
	//flow step size
	double epsilon = 0.02;
	//instantiate the lattice and action
	Lattice lattice(NrDims, shape);
	Wilson action(beta);
	//run through the saved gauge configurations
	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		//load the gauge configuration into U
		SU3_field U(lattice, extdofs);
		std::string ensembleNum = std::to_string(i);
		U.loadSU3FromFile(beta, "Heatbath_4_ORperHB", ensembleNum);
		//instantiate the Gradient flow object
		GradientFlow flowing(action, U, epsilon);
		//instantiate and include the Topological Charge
		TopologicalCharge topCharge(U);
		flowing.Include_TopCharge(topCharge);
		//instantiate and include the Energy density
		//EnergyDensity Edensity(U);
		//flowing.Include_EnergyDensity(Edensity);

		//flow the configuration
		for (int i = 0; i < flowSteps; i++) {
			flowing.flow();
		}
		//save the observables to files
		topCharge.saveTopologicalChargeToFile(beta,flowing.getupdateMethod() ,"Heatbath_4_ORperHB_Flowed" + ensembleNum);
		//Edensity.saveEnergyDensityToFile(beta, ensembleNum + "flowed");
		U.saveSU3ToFile(beta, flowing.getupdateMethod(), "Heatbath_4_ORperHB_Flowed" + ensembleNum);
	}
}

void GenerateLHMCGaugeConfigurations(int NrDims, int extdofs, int shape[], double beta, int ConfigurationStart, int ConfigurationStop, int NrUpdates) {
	
	double LHMCepsilon = 0.33;
	int NrLeaps = 3;

	Lattice lattice(NrDims, shape);
	Wilson action(beta);

	for (int i = ConfigurationStart; i <= ConfigurationStop; i++) {
		std::string ensembleNum = "ensemble_"+ std::to_string(i);
		SU3_field U(lattice, extdofs);
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
		U.saveSU3ToFile(beta, ensembleNum);
		topCharge.saveTopologicalChargeToFile(beta, ensembleNum +"LHMC");
		Edensity.saveEnergyDensityToFile(beta, ensembleNum+"LHMC");
	}
}
void testGradientFlow(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	mpi_debug_breakpoint
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 8,8,8,8 };
	double beta = 6.0;
	double epsilon = 0.01;
	double HMCepsilon = 0.25;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	Wilson action(beta);
	TopologicalCharge topCharge(U);

	LHMC updater(U, action, HMCepsilon);
	for (int i = 0; i < 50; i++) {
		updater.sweep();
		topCharge.calculate(0.0);
	}

	GradientFlow flowing(action, U,epsilon);
	flowing.Include_TopCharge(topCharge);
	topCharge.calculate(0.0);
	for (int i = 0; i < 100; i++) {
		flowing.flow();
	}

	mpiWrapper::end_parallelSession();
}

void testLHMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 8,8,8,8};
	double beta = 6.2;
	double epsilon = 0.5;
	int NrLeaps = 2;
	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	Wilson action(beta);
	mpi_debug_breakpoint
	LHMC updater(U, action, epsilon,NrLeaps);
	//TopologicalCharge topCharge(U);
	//topCharge.calculate(0.0);
	EnergyDensity Edensity(U);
	for (int i = 0; i < 100; i++) {
		updater.sweep();
		//topCharge.calculate(0.0);
		Edensity.calculate(0.0);
	}
	std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";

	mpiWrapper::end_parallelSession();
}
void testSaveLoad(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 2;
	int extdofs = 1;
	int shape[] = { 10,10 };
	Lattice lattice(NrDims, shape);
	SU3_field field(lattice, extdofs);
	su3_mat unit;
	C_double Iunit(0.0, 1.0);
	unit.setToIdentity();
	if (mpiWrapper::id() == 0) {
		for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
			for (int mu = 0; mu < extdofs; mu++) {
				field(i, mu) = i * unit + i * unit * Iunit;
			}
		}
	}
	if (mpiWrapper::id() == 1) {
		for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
			for (int mu = 0; mu < extdofs; mu++) {
				field(i, mu) = 2 * i * unit + 2 * i * unit * Iunit;
			}
		}
	}

	// field.saveSU3ToFile();

	SU3_field field2(lattice, extdofs);
	//field2.loadSU3FromFile();
	if (mpiWrapper::id() == 0) {
		for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
			for (int mu = 0; mu < extdofs; mu++) {
				std::cout << field2(i, mu) << "\n";
			}
		}
	}

	mpiWrapper::end_parallelSession();
}
double getLatticeConstant(double beta) {
	return 0.5 * exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0) * (beta - 6.0) - 0.6667 * (beta - 6.0) * (beta - 6.0) * (beta - 6.0));
}