// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>


int main(int argc, char** argv) {
	//testHMC(argc, argv);
	//testLHMC(argc, argv);
	//testLMC(argc, argv);
	//testSU3_Heatbath(argc, argv);
	//testGradientFlow(argc, argv);
	GenerateGaugeEnsembles(argc, argv);
	//FlowSavedGaugeEnsembles(argc, argv);
	//mpiWrapper::begin_parallelSession(argc, argv);
	//mpi_debug_breakpoint
	//int NrDims = 4;
	//int extdofs = 4;
	//int shape[] = { 8,8,8,8 };
	//double beta = 6.0;
	//double epsilon = 0.01;
	//double HMCepsilon = 0.25;

	//Lattice lattice(NrDims, shape);
	//SU3_field U(lattice, extdofs);
	//U.InitializeHotStart();
	//Wilson action(beta);
	//TopologicalCharge topCharge(U);
	//topCharge.calculate(0.0);
	////std::cout << (U.clover_avg(0, 0,1)* U.clover_avg(0, 0, 1)).ReTr() <<"\n";
	////std::cout << (U.clover_avg1(0, 0,1)* U.clover_avg1(0, 0, 1)).ReTr() << "\n";
	//mpiWrapper::end_parallelSession();


	return 0;
}
void testSU3_Heatbath(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	mpi_debug_breakpoint
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 24,24,24,24 };
	double beta = 0.2;

	Lattice lattice(NrDims, shape);
	//Wilson action(beta);

	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	TopologicalCharge topCharge(U);

	SU3_Heatbath heatbath(U);
	int NrRuns = 2;
	for (int i = 0; i < NrRuns; i++) {
		std::filesystem::path ensembleNum("ensemble_Heatbath");
		heatbath.update(U,30,beta);
		topCharge.calculate(0.0);
		//topCharge.saveTopologicalChargeToFile(6.0, ensembleNum);
		std::filesystem::path ensemblename("ensemble_Heatbath" + std::to_string(i));
		U.saveSU3ToFile(beta,ensemblename);
	}

	mpiWrapper::end_parallelSession();
}
void FlowSavedGaugeEnsembles(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 12,12,12,12 };
	double beta = 6.0;
	double epsilon = 0.01;

	Lattice lattice(NrDims, shape);
	Wilson action(beta);
	int NrEnsembles = 150;
	for (int i = 31; i < NrEnsembles; i++) {
		std::filesystem::path ensembleNum("ensemble_" + std::to_string(i));
		SU3_field U(lattice, extdofs);
		U.loadSU3FromFile(beta, ensembleNum);
		TopologicalCharge topCharge(U);
		EnergyDensity Edensity(U);
		GradientFlow flowing(action, U, epsilon);
		flowing.Include_TopCharge(topCharge);
		flowing.Include_EnergyDensity(Edensity);
		for (int i = 0; i < 600; i++) {
			flowing.flow();
		}
		topCharge.saveTopologicalChargeToFile(beta, ensembleNum);
		Edensity.saveEnergyDensityToFile(beta, ensembleNum);
	}
	//std::filesystem::path specifier("flowed");
	//U.saveSU3ToFile(specifier);
	mpiWrapper::end_parallelSession();
}
void GenerateGaugeEnsembles(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	mpi_debug_breakpoint
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 12,12,12,12 };
	double beta = 6.0;
	
	double LHMCepsilon = 0.33;
	int NrLeaps = 3;

	Lattice lattice(NrDims, shape);
	Wilson action(beta);
	int NrEnsembles = 150;
	for (int i = 45; i < NrEnsembles; i++) {
		std::filesystem::path ensembleNum("ensemble_"+ std::to_string(i));
		SU3_field U(lattice, extdofs);
		U.InitializeHotStart();
		LHMC updater(U, action, LHMCepsilon,NrLeaps);
		for (int i = 0; i < 120; i++) {
			updater.sweep();
		}
		std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";
		MPI_Barrier(mpiWrapper::comm());
		U.saveSU3ToFile(beta, ensembleNum);
	}

	/// //////////////////////////////////////////////////////////////////
	for (int i = 31; i < NrEnsembles; i++) {
		std::filesystem::path ensembleNum("ensemble_" + std::to_string(i));
		SU3_field U(lattice, extdofs);
		U.loadSU3FromFile(beta, ensembleNum);
		TopologicalCharge topCharge(U);
		EnergyDensity Edensity(U);
		GradientFlow flowing(action, U, 0.01);
		flowing.Include_TopCharge(topCharge);
		flowing.Include_EnergyDensity(Edensity);
		for (int i = 0; i < 600; i++) {
			flowing.flow();
		}
		topCharge.saveTopologicalChargeToFile(beta, ensembleNum);
		Edensity.saveEnergyDensityToFile(beta, ensembleNum);
	}


	/// REMOVE
	/// //////////////////////////////////////////////////////////////////

	mpiWrapper::end_parallelSession();
}
void testGradientFlow(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	mpi_debug_breakpoint
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 12,12,12,12 };
	double beta = 6.0;
	double epsilon = 0.01;
	double HMCepsilon = 0.25;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	Wilson action(beta);
	TopologicalCharge topCharge(U);

	LHMC updater(U, action, HMCepsilon);
	for (int i = 0; i < 100; i++) {
		updater.sweep();
		topCharge.calculate(0.0);
	}

	GradientFlow flowing(action, U,epsilon);
	flowing.Include_TopCharge(topCharge);
	topCharge.calculate(0.0);
	for (int i = 0; i < 10000; i++) {
		flowing.flow();
	}

	mpiWrapper::end_parallelSession();
}


void testLMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 4,4,4,4 };
	double beta = 3.1;
	double epsilon = 0.05;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	Wilson action(beta);
	LMC updater(U,epsilon);
	TopologicalCharge topCharge(U);
	std::cout << "updating...\n";
	for (int i = 0; i < 500; i++) {
		updater.update();
		topCharge.calculate(0.0);
	}
	std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";

	mpiWrapper::end_parallelSession();
}



void testHMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	mpi_debug_breakpoint
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 6,6,6,6};
	double beta = 6.0;
	double epsilon = 0.007;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	U.InitializeHotStart();
	Wilson action(beta);
	HMC updater(U, action, epsilon);
	TopologicalCharge topCharge(U);
	topCharge.calculate(0.0);
	for (int i = 0; i < 150; i++) {
		updater.update();
		topCharge.calculate(0.0);
	}
	std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";

	mpiWrapper::end_parallelSession();
}

void testLHMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);
	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 12,12,12,12};
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
