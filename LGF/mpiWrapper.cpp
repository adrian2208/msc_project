#include "mpiWrapper.h"

int mpiWrapper::proc_id = 0;
int mpiWrapper::Nprocs = 1;
MPI_Comm mpiWrapper::Comm;

void mpiWrapper::begin_parallelSession(int argc, char** argv) {
	Comm = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(Comm, &(Nprocs));
	MPI_Comm_rank(Comm, &(proc_id));
	
	if (proc_id == 0) {
		std::cout << "Started Parallel Session with " << Nprocs << " processes...\n";
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

MPI_Comm mpiWrapper::comm(){
	return Comm;
}

void mpiWrapper::mpi_openFile(MPI_File &file, const char* filename) {
	MPI_File_open(Comm, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
}
void mpiWrapper::mpi_closeFile(MPI_File& file) {
	MPI_File_close(&file);
}