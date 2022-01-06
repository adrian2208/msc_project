// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
#include <filesystem>
void print_array(int** array, int nRows, int nCols);
void print_array(int** array, int nRows, int nCols) {
	for (int i = 0; i < nRows; ++i) {
		for (int j = 0; j < nCols; ++j) {
			std::cout << array[i][j] << ' ';
		}
		std::cout << std::endl;
	}
}
void print_vec(int* array, int nitems);
void print_vec(int* array, int nitems) {
	for (int j = 0; j < nitems; ++j) {
			std::cout << array[j] << ' ';
	}
	std::cout << std::endl;
}


int main(int argc, char **argv){
	
	mpiWrapper::begin_parallelSession(argc, argv);

	int NrDims = 2;
	int extdofs = 1;
	int shape[] = {100,100};
	Lattice lattice(NrDims, shape);
	SU3_field field(lattice, extdofs);
	su3_mat unit;
	unit.setToIdentity();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < extdofs; mu++) {
			field(i, mu) = i*unit;
		}
	}
	field.saveSU3ToFile();

	SU3_field field2(lattice, extdofs);
	field2.loadSU3FromFile();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < extdofs; mu++) {
			std::cout << field2(i, mu) << "\n";
		}
	}

	mpiWrapper::end_parallelSession();
	
	return 0;
}
