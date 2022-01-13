// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>
#define DEBUG

int main(int argc, char** argv) {
	//mpiWrapper::begin_parallelSession(argc, argv);
	//double epsilon = 0.0001;
	//Random rng;
	//su3_mat mat;
	//SU3_gen generators;
	//int number_accepted = 0;
	//for (int i = 0; i < 10; i++) {
	//	//rng.rnd_su3_alg(mat,epsilon);
	//	mat.setToZeros();
	//	for (int j = 0; j < 8; j++) {
	//		mat = mat + rng.Uniform_Double() * generators(j);
	//	}
	//	mat = 100*mat.timesMinusI();
	//	number_accepted += IsHermTrLess(mat);
	//	HermTrLessExp(mat);
	//	std::cout << mat << "\n";
	//}
	//std::cout << "number accepted: " << number_accepted << "\n";
	//mpiWrapper::end_parallelSession();
	testHMC(argc, argv);

	return 0;
}
void testHMC(int argc, char** argv) {
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 4;
	int extdofs = 4;
	int shape[] = { 2,2,2,2 };
	double beta = 1.0;
	double epsilon = 0.1;

	Lattice lattice(NrDims, shape);
	SU3_field U(lattice, extdofs);
	su3_mat unit;
	unit.setToIdentity();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < extdofs; mu++) {
			U(i, mu) = unit;
		}
	}
	Wilson action(beta);
	HMC updater(U, action, epsilon);
	updater.update();

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
