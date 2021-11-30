// LGF.cpp : Defines the entry point for the application.
//C:\Users\adria\Documents\msc_project\LGF\out\build\x64-Release>mpiexec.exe -n 6 LGF.exe

#include "LGF.h"



int main(int argc, char **argv){

	mpiWrapper::begin_parallelSession(argc, argv);
	int shape[] = { 4,4 };
	Lattice lattice(2, shape);
	mpiWrapper::end_parallelSession();
	
	return 0;
}
