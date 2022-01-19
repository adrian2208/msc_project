﻿// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>

int main(int argc, char** argv) {
	testHMC(argc, argv);
	//testLMC(argc, argv);

	//mpiWrapper::begin_parallelSession(argc, argv);
	//int NrDims = 4;
	//int extdofs = 4;
	//int shape[] = { 2,2,2,2 };
	//double beta = 6.1;
	//double epsilon = 0.01;

	//Wilson action(beta);
	//Lattice lattice(NrDims, shape);
	//SU3_field U(lattice, extdofs);
	//SU3_field F(lattice, extdofs);
	//su3_mat unit;
	//Random rng;
	//su3_mat LinCombGen;
	//SU3_gen generators;
	//unit.setToIdentity();
	//for (int attempts = 0; attempts < 50; attempts++) {
	//	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
	//		for (int mu = 0; mu < extdofs; mu++) {
	//			//U(i, mu) = unit;
	//			for (int j = 0; j < 8; j++) {
	//				LinCombGen = LinCombGen + rng.Gaussian_Double(0.0, 1.0) * generators(j);
	//			}
	//			LinCombGen = LinCombGen.timesI();
	//			U(i, mu) = HermTrLessExp(LinCombGen);
	//			isSpecialUnitary(U(i, mu), false);
	//			LinCombGen.setToZeros();
	//			
	//		}
	//	}
	//	action.calculate_Force(U,F);
	//	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
	//		for (int mu = 0; mu < extdofs; mu++) {
	//			std::cout << F(i, mu) << "\n\n";
	//		}
	//		std::cout << "\n";
	//	}

	//	std::cout << "avg. plaquette = " << U.Avg_Plaquette() << "\n";
	//	std::cout << "Action = " << action.calculate_Action(U) << "\n";
	//	std::cout << "normalized Action = " << action.calculate_Action(U)/((double)U.getLatticePtr().m_responsible_Volume*12.0) << "\n\n";
	//}
	//mpiWrapper::end_parallelSession();


	return 0;
}
void testLMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 6,6,6,6 };
	double beta = 6.1;
	double epsilon = 0.05;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	su3_mat unit;
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	unit.setToIdentity();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < extdofs; mu++) {
			U(i, mu) = unit;
			//for (int j = 0; j < 8; j++) {
			//	LinCombGen = LinCombGen + rng.Gaussian_Double(0.0, 1.0) * generators(j);
			//}
			//LinCombGen = LinCombGen.timesI();
			//U(i, mu) = HermTrLessExp(LinCombGen);
			//LinCombGen.setToZeros();
		}
	}

	Wilson action(beta);
	LMC updater(U,epsilon);
	std::cout << "updating...\n";
	for (int i = 0; i < 100; i++) {
		updater.update();
	}
	std::cout << "acceptance rate1: " << updater.acceptanceRate() << "\n";

	mpiWrapper::end_parallelSession();
}


void testHMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 4,4,4,4};
	double beta = 3.0;
	double epsilon = 0.01;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	su3_mat unit;
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	unit.setToIdentity();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < extdofs; mu++) {
			//U(i, mu) = unit;
			for (int j = 0; j < 8; j++) {
				LinCombGen = LinCombGen + rng.Gaussian_Double(0.0, 1.0) * generators(j);
			}
			LinCombGen = LinCombGen.timesI();
			U(i, mu) = HermTrLessExp(LinCombGen);
			LinCombGen.setToZeros();
		}
	}

	Wilson action(beta);
	HMC updater(U, action, epsilon);
	for (int i = 0; i < 10; i++) {
		updater.update();
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

	field.saveSU3ToFile();

	SU3_field field2(lattice, extdofs);
	field2.loadSU3FromFile();
	if (mpiWrapper::id() == 0) {
		for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
			for (int mu = 0; mu < extdofs; mu++) {
				std::cout << field2(i, mu) << "\n";
			}
		}
	}

	mpiWrapper::end_parallelSession();
}
