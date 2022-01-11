// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>

int main(int argc, char** argv) {
	//testSaveLoad(argc, argv);
	//su3_mat unit;
	//unit.setToIdentity();
	//unit[8] = C_double(-2.0, 0.0);
	//unit = unit * (1.0 / sqrt(3));
	//unit = unit.timesI();
	//HermTrLessExp(unit);
	//std::cout << unit << "\n";
	Random rand;
	su3_mat unit;
	rand.rnd_su3_alg(unit, 1.0);
	//unit[8] = C_double(-2.0, 0.0);
	//unit = unit * (1.0 / sqrt(3));
	//std::cout << unit << "\n\n\n";
	//std::cout << unit.at() << "\n";
	//std::cout << "trace: " << unit.Tr() << "\n";
	//std::cout << "H+H.dagger: " << unit+unit.dagger() << "\n";
	return 0;
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
