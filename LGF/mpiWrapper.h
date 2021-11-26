#pragma once
#include "mpi.h"

class mpiWrapper {
private:
	int argc;
	char** argv;
	int proc_id;
	int Nprocs;
public:
	void begin_parallelSession(int argc, char** argv) {
		MPI_Comm Comm= MPI_COMM_WORLD;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(Comm, &(proc_id));
		MPI_Comm_size(Comm, &(Nprocs));
	}
	
};
mpiWrapper com;