#include "Lattice.h"

Lattice::Lattice(int Ndims, int shape[]) {
	
	//Storing local variables
	m_Ndims = Ndims;
	m_shape = new int[Ndims];
	m_coordinate = new int[Ndims];
	m_coorFwd = new int[Ndims];
	m_coorBack = new int[Ndims];
	//Calculates the total volume of the lattice & writing to member variables
	m_totalVolume = 1;
	m_thisProc_Volume = 0;
	m_responsible_Volume = 0;
	for (int i = 0; i < Ndims; i++) {
		m_shape[i] = shape[i];
		m_totalVolume *= shape[i];
		m_coordinate[i] = 0;
	}

	partition_lattice();

}

//iterates through an N- dimensional array of arbitrary shape
inline void traverse_lattice(int* coordinate, int nDims, int* shape) {
	for (int i = nDims - 1; i >= 0; i--) {							
		if (coordinate[i] < shape[i]-1) { 
			coordinate[i]++;	
			break;
		}									
		else {								
			coordinate[i] = 0;			
		}									
	}
}
//assigns even (0), or odd (1) parity to a coordinate
inline int parity(int* coordinate,int nDims) {
	int parity = 0;
	for (int i = 0; i < nDims; i++) {
		parity += coordinate[i];
	}
	return (parity % 2);
}

void Lattice::partition_lattice(){
	m_thisProc_TotalIndex = new int[m_totalVolume];
	int totalIdx;
	int id;
	
	int i=0;
	while (i < m_totalVolume) {//works with traverse_lattice, iterating through the lattice
																//__
		totalIdx = totalIndex(m_coordinate);					//	|
		id = Coordinate_ProcID(m_coordinate);					//  |
		if (mpiWrapper::id() == id) {							//  |-- identifies whether the working process
			m_thisProc_TotalIndex[m_thisProc_Volume] = totalIdx;//  |	is responsible for that coordinate, 
			m_thisProc_Volume++;								//  |	incrementing the proc volumes if it is.
			m_responsible_Volume++;								//	|
		}														//__|
																		//__
		else if(is_SharedMemory(m_coordinate, m_coorFwd, m_coorBack)) {	//	|-- If the process is not responsible,
			m_thisProc_TotalIndex[m_thisProc_Volume] = totalIdx;		//	|	but that coordinate is within the
			m_thisProc_Volume++;										//	|	neighbouring process' shared memory
		}																//__|	region, increment the proc volume anyway.
		//returns new m_coordinate by reference eventually traversing the lattice															
		traverse_lattice(m_coordinate, m_Ndims, m_shape); 												
		i++;
	}

	//Initializing arrays for the starting and stopping internal indices 
	//according to the ID of the process responsible and the parity of the site.
	//(Note that not every process needs to have every other process as a neighbour, ...
	//this solution is just simpler than counting the number of unique neighbouring processes)
	m_InternalIdx_start = new int* [mpiWrapper::nProcs()];
	m_InternalIdx_stop = new int* [mpiWrapper::nProcs()];
	for (i = 0; i < mpiWrapper::nProcs(); i++) {
		m_InternalIdx_start[i] = new int[2];
		m_InternalIdx_stop[i] = new int[2];
	}
	m_InternalIdx_start[0][0] = 0;
	m_InternalIdx_stop[0][0] = 0;
	//Initializing array for converting between the total index across all processes
	//and processes' internal lattice indices
	m_TotalToInternal_idx = new int[m_totalVolume];
	for (i = 0; i < m_totalVolume; i++) {
		m_TotalToInternal_idx[i] = INT_MAX;
	}
	//Initializing array for converting between processes' internal lattice indices
	//and the total index across all processes
	m_InternalToTotal_idx = new int[m_thisProc_Volume];
	//Initializing array for converting between processes' internal lattice indices
	//and the ID of the process responsible
	m_Internal_IdxToProcId = new int[m_thisProc_Volume];
	//Initializing array for converting between processes' internal lattice indices
	//and the parity of that site
	m_Internal_IdxToParity = new int[m_thisProc_Volume];
	int InternalIdx;
	//The internal index for the sites stored by this process:
	//1: prioritizes the smallest process number
	//2: prioritizes even parity over odd parity.
	for (int proc_id = 0; proc_id < mpiWrapper::nProcs(); proc_id++) {
		//avoid the one case that would make m_InternalIdx_stop[proc_id - 1][1] an illegal memory read
		if (proc_id != 0) {
			//the start index of even parity sites of the next process
			//begins at the stop index of odd parity sites of the former process
			m_InternalIdx_start[proc_id][0] =  m_InternalIdx_stop[proc_id - 1][1];
			//beginning from the start index, the stop index will later be incremented
			m_InternalIdx_stop[proc_id][0] = m_InternalIdx_start[proc_id][0];
		}
		for (int par = 0; par <= 1; par++) {
			//execute only after m_InternalIdx_stop[proc_id][0] has aquired the correct value
			if (par == 1) {
				//the start index of odd parity sites of the same process
				//begins at the stop index of even parity sites of the same process
				m_InternalIdx_start[proc_id][1] = m_InternalIdx_stop[proc_id][0];
				//beginning from the start index, the stop index will later be incremented
				m_InternalIdx_stop[proc_id][1] = m_InternalIdx_start[proc_id][1];
			}
			for (int j = 0; j < m_thisProc_Volume; j++) {
				//update the working total coordinate to the next point stored by the process
				IndexToCoordinate(m_thisProc_TotalIndex[j], m_coordinate);
				if (proc_id == Coordinate_ProcID(m_coordinate) && par == parity(m_coordinate, m_Ndims)) {
					InternalIdx = m_InternalIdx_stop[proc_id][par];
					//Which process is responsible for the point of that internal index
					m_Internal_IdxToProcId[InternalIdx] = proc_id;
					//What is the parity of the point with that internal index
					m_Internal_IdxToParity[InternalIdx] = par;
					//provides maps between the total and internal indices
					m_InternalToTotal_idx[InternalIdx] = m_thisProc_TotalIndex[j];
					m_TotalToInternal_idx[m_thisProc_TotalIndex[j]] = InternalIdx;
					
					//Increment the stop index for that process and parity if a site
					//of that process and parity is found
					m_InternalIdx_stop[proc_id][par]++;
				}
	
			}
		}
	}

	//Initializing arrays for storing nearest neighbour information
	m_cntr = new int* [m_thisProc_Volume];
	m_fwd = new int* [m_thisProc_Volume];
	m_back = new int* [m_thisProc_Volume];
	for (i = 0; i < m_thisProc_Volume; i++) {
		m_cntr[i] = new int[m_Ndims];
		m_fwd[i] = new int[m_Ndims];
		m_back[i] = new int[m_Ndims];
	}
	for (i = 0; i < m_thisProc_Volume; i++) {
		IndexToCoordinate(m_InternalToTotal_idx[i], m_coordinate);
		for (int j = 0; j < m_Ndims; j++) {
			m_cntr[i][j] = m_coordinate[j];
			NearestNeighbour(m_coordinate, j, m_coorFwd, m_coorBack);
			//if (m_Internal_IdxToProcId[i] == mpiWrapper::id()) {
				m_fwd[i][j] = m_TotalToInternal_idx[totalIndex(m_coorFwd)];
				m_back[i][j] = m_TotalToInternal_idx[totalIndex(m_coorBack)];
			//}
			//else{
			//	if (m_TotalToInternal_idx[totalIndex(m_coorFwd)] != INT_MAX) {
					
			//	}
		}
	}

	DistributeBuffers();
	

}

void Lattice::DistributeBuffers() {
	m_BufferSize = new int [mpiWrapper::nProcs()];
	m_Buffer_receive = new int* [mpiWrapper::nProcs()];

	MPI_Request request;
	MPI_Status status;
	int SizeofProcDomain;
	int* TotalIdx_procDomain;
	int SendingTo_id;//The ID of the process which will be recieving information
	int RecFrom_id;//The ID of the process from which information will be recieved
	
	for (int i = 1; i < mpiWrapper::nProcs(); i++) {
		//recieving and sending from two different processes. Never equal to id.
		SendingTo_id = (mpiWrapper::id() + i) % mpiWrapper::nProcs();
		RecFrom_id = (mpiWrapper::id() + mpiWrapper::nProcs() - i) % mpiWrapper::nProcs();
		//The size of the package <<SendingTo_id>> should expect to receive based on the amount of sites
		//shared between <<ID>> and <<SendingTo_id>> 
		SizeofProcDomain = m_InternalIdx_stop[SendingTo_id][1] - m_InternalIdx_start[SendingTo_id][0];

		//Communicate - wait for resolution
		MPI_Isend(&SizeofProcDomain, 1, MPI_INT, SendingTo_id, mpiWrapper::id() * mpiWrapper::nProcs() + SendingTo_id, mpiWrapper::comm(), &request);
		MPI_Recv(&m_BufferSize[RecFrom_id], 1, MPI_INT, RecFrom_id, RecFrom_id * mpiWrapper::nProcs() + mpiWrapper::id(), mpiWrapper::comm(), &status);
		MPI_Wait(&request, &status);

		//for each shared site belonging to process <<SendingTo_id>>, convert the internal parameterization (index)
		//to the total one, which is shared by all processes
		TotalIdx_procDomain = new int[SizeofProcDomain];
		for (int j = 0; j < SizeofProcDomain; j++) {
			TotalIdx_procDomain[j] = m_InternalToTotal_idx[m_InternalIdx_start[SendingTo_id][0] + j];
			//std::cout << TotalIdx_procDomain[j] << " ";
			//std::cout << TotalIdx_procDomain[j] << " ";
		}
		//std::cout << "\n";

		//send the list of shared sites belonging to the process in their total index to the process in question
		MPI_Isend(TotalIdx_procDomain, SizeofProcDomain, MPI_INT, SendingTo_id, mpiWrapper::id() * mpiWrapper::nProcs() + SendingTo_id, mpiWrapper::comm(), &request);
		//prepare the array which will contain the list of ID's sites of responsibility contained in the domain of <<RecFrom_id>>
		m_Buffer_receive[RecFrom_id] = new int[m_BufferSize[RecFrom_id]];
		MPI_Recv(m_Buffer_receive[RecFrom_id], m_BufferSize[RecFrom_id], MPI_INT, RecFrom_id, RecFrom_id * mpiWrapper::nProcs() + mpiWrapper::id(), mpiWrapper::comm(), &status);
		//The sites are currently given by their total index and have to be assigned their internal parameterization
		for (int j = 0; j < m_BufferSize[RecFrom_id]; j++) {
			
			m_Buffer_receive[RecFrom_id][j] = m_TotalToInternal_idx[m_Buffer_receive[RecFrom_id][j]];
			//std::cout << m_Buffer_receive[RecFrom_id][j] << " ";
		}
		//std::cout << "From: " << RecFrom_id << "\n";
		MPI_Wait(&request, &status);
		delete[] TotalIdx_procDomain;
	}
}


int Lattice::Coordinate_ProcID(int* coordinate){
	int NrProcs = mpiWrapper::nProcs();
	int remainder = (int)m_shape[0] % NrProcs;
	if (remainder == 0 && NrProcs>1) {
		return coordinate[0] / (m_shape[0] / NrProcs);
	}
	return 0;
}


int Lattice::totalIndex(int *coordinate){
	//  *(0,1)--*(1,1)
	//	|		|		: Coordinate
	//	*(0,0)--*(1,0)
	//
	//	1--3
	//	|  |			: totalIndex
	//	0--2
	int idx = 0;
	for (int i = 0; i < m_Ndims - 1; i++) {
		idx = m_shape[i + 1] * (idx + coordinate[i]);
	}
	return idx+coordinate[m_Ndims-1];
}

void Lattice::IndexToCoordinate(int index, int* coordinate) {
	for (int i = m_Ndims - 1; i > 0; i--) {
		coordinate[i] = index % m_shape[i];
		index = (index - coordinate[i]) / m_shape[i];
	}
	coordinate[0] = index;
}

// This virtual method encodes periodic boundary conditions for all directions on the lattice
void Lattice::NearestNeighbour(int* coordinate, int direction, int* fwd_coor, int* back_coor){
	for (int i = 0; i < m_Ndims; i++) {
		if (i == direction) {
			fwd_coor[i] = (coordinate[i] + 1) % m_shape[i];
			back_coor[i] = (coordinate[i] - 1 + m_shape[i]) % m_shape[i];
		}
		else {
			fwd_coor[i] = coordinate[i];
			back_coor[i] = coordinate[i];
		}
	}
}

bool Lattice::is_SharedMemory(int* coordinate, int* coorFwd, int* coorBack){
	int id_fwd;
	int id_back;
	for (int i = 0; i < m_Ndims; i++) {
		NearestNeighbour(coordinate, i, coorFwd, coorBack);
		id_fwd = Coordinate_ProcID(coorFwd);
		id_back = Coordinate_ProcID(coorBack);
		if (mpiWrapper::id() == id_fwd || mpiWrapper::id() == id_back) {
			return true;
		}
	}
	return false;
}

void Lattice::print_indices(){
	for (int i = 0; i < m_thisProc_Volume; i++) {
			std::cout << m_thisProc_TotalIndex[i] << ",";
	}
	std::cout << "\n";
	MPI_Barrier(mpiWrapper::comm());
}

void Lattice::print() {
	int** myArray = new int* [m_shape[0]];
	for (int i = 0; i < m_shape[0]; i++) {
		myArray[i] = new int[m_shape[1]];
	}

	int* coor = new int[2];
	coor[0] = 0;
	coor[1] = 0;

	int i = 0;
	while (i < m_totalVolume) {
		int totalIdx = totalIndex(coor);
		int x = coor[0];
		int y = coor[1]; 
		myArray[x][y]=totalIdx;
		std::cout << x << "," << y <<"," << totalIdx << "\n";
		traverse_lattice(coor, m_Ndims, m_shape);
		
		i++;
	}

	for (i = 0; i < m_shape[0]; ++i){
		for (int j = 0; j < m_shape[1]; ++j){
			std::cout << myArray[i][j] << ' ';
		}
		std::cout << std::endl;
	}

	MPI_Barrier(mpiWrapper::comm());
}
