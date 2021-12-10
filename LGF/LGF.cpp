// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"
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
	int shape[] = { 4,4};
	Lattice lattice(2, shape);
	Field<int> field(lattice, 1);
	mpiWrapper::end_parallelSession();
	

	
	return 0;
}
