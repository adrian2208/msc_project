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
	int shape[] = { 64,64,64,64};
	Lattice lattice(4, shape);
	SU3_field field(lattice, 4);
	//field.loadFromFile();
	su3_mat unit;
	unit.setToIdentity();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < 1; mu++) {
			field(i, mu) = unit;
		}
	}
	field.saveSU3ToFile();
	field.loadSU3FromFile();
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < 1; mu++) {
			std::cout << field(i, mu) << "\n";
		}
	}
	mpiWrapper::end_parallelSession();
	
	/*Random rand;
	su3_mat mat;
	double epsilon = 1;
	rand.rnd_su3_alg(mat,epsilon);

	std::cout << mat.dagger().det() << "\n";
	std::cout << mat.det() << "\n";
	std::cout << mat* mat.dagger() << "\n";
	*/
	return 0;
}
