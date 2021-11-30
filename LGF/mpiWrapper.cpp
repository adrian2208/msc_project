#include "mpiWrapper.h"

int mpiWrapper::proc_id = 0;
int mpiWrapper::Nprocs = 1;
MPI_Comm mpiWrapper::Comm;

void mpiWrapper::begin_parallelSession(int argc, char** argv) {
	Comm = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(Comm, &(Nprocs));
	MPI_Comm_rank(Comm, &(proc_id));
	
	for (int i = 0; i < Nprocs; i++) {
		if (i == proc_id) {
			std::cout << "Process ID: " << i << "\n" << "Nr. Processes: " << Nprocs << "\n";
			MPI_Barrier(Comm);
		}
	}
}
void mpiWrapper::end_parallelSession() {
	MPI_Finalize();
}
int mpiWrapper::id() {
	return proc_id;
}
int mpiWrapper::nProcs() {
	return Nprocs;
}