#pragma once
#include "mpi.h"
#include <iostream>

#ifndef MPIWRAPPER_INCLUDED
#define MPIWRAPPER_INCLUDED
class mpiWrapper {
private:
	static MPI_Comm Comm;
	int argc;
	char** argv;
	static int proc_id;
	static int Nprocs;

public:

	static void begin_parallelSession(int argc, char** argv);
	static void end_parallelSession();
	static int id();
	static int nProcs();
	static MPI_Comm comm();

	static void mpi_openFile(MPI_File& file, const char* filename);
	static void mpi_closeFile(MPI_File& file);

};

#endif