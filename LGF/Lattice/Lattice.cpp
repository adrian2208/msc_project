#include "Lattice.h"
#include <assert.h>
#include <climits>
Lattice::Lattice(SimulationParameters& params,bool cloversRequired) {
	
	//Storing local variables
	m_cloversRequired = cloversRequired;
	type = params.getLatticeType();
	m_Ndims = params.getNrDims();
	m_shape = new int[m_Ndims];
	m_cuts = new int[m_Ndims];
	m_coordinate = new int[m_Ndims];
	m_coorFwd = new int[m_Ndims];
	m_coorBack = new int[m_Ndims];

	m_coorFwdFwd = new int[m_Ndims];
	m_coorBackBack = new int[m_Ndims];
	m_coorFwdBack = new int[m_Ndims];
	m_coorBackFwd = new int[m_Ndims];

	m_coorBackFwdFwd = new int[m_Ndims];//
	m_coorBackBackFwd = new int[m_Ndims];//
	m_coorFwdFwdFwd = new int[m_Ndims];//
	m_coorBackBackBack = new int[m_Ndims];//
	m_coorBackFwdBack = new int[m_Ndims];//
	m_coorFwdFwdBack = new int[m_Ndims];//
	

	m_coorFwdFwdFwdFwd = new int[m_Ndims];//NOT IN USE
	m_coorBackBackBackBack = new int[m_Ndims];//
	m_coorFwdFwdFwdBack = new int[m_Ndims];//
	m_coorBackBackBackFwd = new int[m_Ndims];//
	m_coorBackFwdBackFwd = new int[m_Ndims];//
	m_coorBackFwdBackBack = new int[m_Ndims];//NOT IN USE
	//Calculates the total volume of the lattice & writing to member variables
	m_totalVolume = 1;
	m_thisProc_Volume = 0;
	m_responsible_Volume = 0;
	for (int i = 0; i < m_Ndims; i++) {
		m_shape[i] = params.getshape(i);
		m_cuts[i] = params.getcuts(i);
		if (m_shape[i] % 2 != 0) {
			if (mpiWrapper::id() == 0) {
				std::cout << "Check that shape matches the number of dimensions!\n";
				std::cout << "Even- Odd scheme requires an even number of lattice points in each direction with the torus topology!\n";
			}
			mpiWrapper::end_parallelSession();
			exit(EXIT_FAILURE);
		}
		m_totalVolume *= m_shape[i];
		m_coordinate[i] = 0;
	}
	//CheckPartitioning();
	//Print_Partitioning();
	
	partition_lattice();
	//Print_Partitioning();
	
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
			//if (mpiWrapper::id() == 0) {
			//	std::cout << "coor :" << m_coordinate[0] << " " << m_coordinate[1] << " " << m_coordinate[2] << " " << m_coordinate[3] << std::endl;
			//}
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
	//m_cntr = new int* [m_thisProc_Volume+1];
	//m_fwd = new int* [m_thisProc_Volume+1];
	//m_back = new int* [m_thisProc_Volume+1];
	//for (i = 0; i < m_thisProc_Volume+1; i++) {
	//	m_cntr[i] = new int[m_Ndims];
	//	m_fwd[i] = new int[m_Ndims];
	//	m_back[i] = new int[m_Ndims];
	//}
	m_cntr = new int* [m_thisProc_Volume];
	m_fwd = new int* [m_thisProc_Volume];
	m_back = new int* [m_thisProc_Volume];

	for (i = 0; i < m_thisProc_Volume; i++) {
		m_cntr[i] = new int[m_Ndims];
		m_fwd[i] = new int[m_Ndims];
		m_back[i] = new int[m_Ndims];
	}
	std::cout << "total Volume: " << m_totalVolume << std::endl;
	for (i = 0; i < m_thisProc_Volume; i++) {
		IndexToCoordinate(m_InternalToTotal_idx[i], m_coordinate);
		for (int j = 0; j < m_Ndims; j++) {
			m_cntr[i][j] = m_coordinate[j];
			NearestNeighbour(m_coordinate, j, m_coorFwd, m_coorBack);
			int idxfwd, idxbck;
			idxfwd = totalIndex(m_coorFwd);
			idxbck = totalIndex(m_coorBack);

			//std::cout << "coor :" << m_coordinate[0] << " " << m_coordinate[1] << " " << m_coordinate[2] << " " << m_coordinate[3] << std::endl;
			//std::cout << "back_coor :" << m_coorBack[0] << " " << m_coorBack[1] << " " << m_coorBack[2] << " " << m_coorBack[3] << std::endl;
			//std::cout << "idxbck: " << idxbck << std::endl;
			//std::cout << "coor :" << m_coordinate[0] << " " << m_coordinate[1] << " " << m_coordinate[2] << " " << m_coordinate[3] << std::endl;
			//std::cout << "fws_coor :" << m_coorFwd[0] << " " << m_coorFwd[1] << " " << m_coorFwd[2] << " " << m_coorFwd[3] << std::endl;
			//std::cout << "idxfwd: " << idxfwd << std::endl;


			if (idxfwd < m_totalVolume) {
				m_fwd[i][j] = m_TotalToInternal_idx[idxfwd];
			}
			else {
				m_fwd[i][j] = m_thisProc_Volume;
			}
			if (idxbck < m_totalVolume) {
				m_back[i][j] = m_TotalToInternal_idx[idxbck];
			}
			else {
				m_back[i][j] = m_thisProc_Volume;
			}
		}
	}
	//for (int j = 0; j < m_Ndims; j++) {
	//	m_fwd[m_thisProc_Volume][j] = m_thisProc_Volume;
	//	m_back[m_thisProc_Volume][j] = m_thisProc_Volume;
	//}
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
		//m_BufferSize[SendingTo_id] = SizeofProcDomain;
		//Communicate - wait for resolution
		MPI_Isend(&SizeofProcDomain, 1, MPI_INT, SendingTo_id, mpiWrapper::id() * mpiWrapper::nProcs() + SendingTo_id, mpiWrapper::comm(), &request);
		MPI_Recv(&m_BufferSize[RecFrom_id], 1, MPI_INT, RecFrom_id, RecFrom_id * mpiWrapper::nProcs() + mpiWrapper::id(), mpiWrapper::comm(), &status);
		MPI_Wait(&request, &status);
		
		//for each shared site belonging to process <<SendingTo_id>>, convert the internal parameterization (index)
		//to the total one, which is shared by all processes
		TotalIdx_procDomain = new int[SizeofProcDomain];
		for (int j = 0; j < SizeofProcDomain; j++) {
			TotalIdx_procDomain[j] = m_InternalToTotal_idx[m_InternalIdx_start[SendingTo_id][0] + j];

		}

		//send the list of shared sites belonging to the process in their total index to the process in question
		MPI_Isend(TotalIdx_procDomain, SizeofProcDomain, MPI_INT, SendingTo_id, mpiWrapper::id() * mpiWrapper::nProcs() + SendingTo_id, mpiWrapper::comm(), &request);
		//prepare the array which will contain the list of ID's sites of responsibility contained in the domain of <<RecFrom_id>>
		m_Buffer_receive[RecFrom_id] = new int[m_BufferSize[RecFrom_id]];
		MPI_Recv(m_Buffer_receive[RecFrom_id], m_BufferSize[RecFrom_id], MPI_INT, RecFrom_id, RecFrom_id * mpiWrapper::nProcs() + mpiWrapper::id(), mpiWrapper::comm(), &status);
		//The sites are currently given by their total index and have to be assigned their internal parameterization
		for (int j = 0; j < m_BufferSize[RecFrom_id]; j++) {
			m_Buffer_receive[RecFrom_id][j] = m_TotalToInternal_idx[m_Buffer_receive[RecFrom_id][j]];
		}
		MPI_Wait(&request, &status);
		delete[] TotalIdx_procDomain;
	}
}

void Lattice::CheckPartitioning() {
	//Check if partitioning is of type TPancake
	TPancake_Flag = true;
	int temp = m_cuts[0];
	for (int i = 1; i < m_Ndims; i++) {
		if (m_cuts[i] != 0) {
			TPancake_Flag = false;
		}
	}


//	STFF_Flag = false;
//	Trivial_Flag = false;
//	
//	int NrProcs_MAX = 1024;
//	int NrProcs = mpiWrapper::nProcs();
//	if (NrProcs == 1) {
//		if (mpiWrapper::id() == 0) {
//			std::cout << "Single process -> Trivial Partitioning" << std::endl;
//		}
//		Trivial_Flag = true;
//		return;
//	}
//	if (mpiWrapper::nProcs() <= m_shape[0] && m_shape[0] % mpiWrapper::nProcs() == 0) {
//		if (mpiWrapper::id() == 0) std::cout << "Processes evenly divide the first dimension -> Temporal Pancake partitioning\n";
//		TPancake_Flag = true;
//		return; 
//	}
//	int NumberOfPartitions = 1;
//	int Divisions = 0;
//	int xDivisions;
//	int yDivisions;
//	
//	if (mpiWrapper::id() == 0) {
//		std::cout << "The following are numbers of processes suitable for this lattice upto " << NrProcs_MAX << std::endl;
//	}
//	while (NumberOfPartitions <= NrProcs_MAX) {
//		NumberOfPartitions = (int)floor((double)((Divisions + 2) * (Divisions + 2)) / 4.0);
//
//		if (Divisions % 2 != 0) {
//			yDivisions = (Divisions - 1) / 2;
//			xDivisions = (Divisions - 1) / 2 + 1;
//		}
//		else {
//			yDivisions = Divisions / 2;
//			xDivisions = Divisions / 2;
//		}
//		if (m_shape[0] % (xDivisions + 1) != 0) {
//			Divisions++;
//			continue;
//		}
//		if (m_shape[1] % (yDivisions + 1) != 0) {
//			Divisions++;
//			continue;
//		}
//		if (mpiWrapper::id() == 0) { 
//			std::cout << NumberOfPartitions << " ";
//		}
//		if (mpiWrapper::nProcs() == NumberOfPartitions) {
//			STFF_Flag = true;
//		}
//		Divisions++;
//	}
//	if (mpiWrapper::id() == 0) {
//		std::cout << std::endl;
//	}
//	if (STFF_Flag) {
//		if (mpiWrapper::id() == 0) {
//			std::cout << "Success! "<< NrProcs << " processes utilized in Space-time French Fries partitioning." << std::endl;
//		}
//		return;
//	}
//	if (mpiWrapper::id() == 0) {
//		std::cout << "Failure! " << NrProcs << " processes could not be distributed on the lattice" << std::endl;
//	}
//	mpiWrapper::end_parallelSession();
//	exit(EXIT_FAILURE);
}


int Lattice::Coordinate_ProcID(int* coordinate){
	
	//if (Trivial_Flag) {
	//	return 0;
	//}
	//if (TPancake_Flag) {
		//int NrProcs = mpiWrapper::nProcs();
		//assert(m_shape[0] % NrProcs == 0);
		//int remainder = m_shape[0] % NrProcs;
		//if (remainder == 0 && NrProcs > 1) {
		//	return coordinate[0] / (m_shape[0] / NrProcs);
		//}
	//}
	//if (STFF_Flag) {
	//	int NrProcs = mpiWrapper::nProcs();
	//	int NumberOfPartitions = 1;
	//	int Divisions = 0;
	//	int xDivisions;
	//	int yDivisions;
	//	while (NumberOfPartitions <= NrProcs) {
	//		NumberOfPartitions = (int)floor((double)((Divisions + 2) * (Divisions + 2)) / 4.0);
	//		Divisions++;
	//	}
	//	Divisions -= 2;
	//	NumberOfPartitions = (int)floor((double)((Divisions + 2) * (Divisions + 2)) / 4.0);
	//	if (Divisions % 2 != 0) {
	//		yDivisions = (Divisions - 1) / 2;
	//		xDivisions = (Divisions - 1) / 2 + 1;
	//	}
	//	else {
	//		yDivisions = Divisions / 2;
	//		xDivisions = Divisions / 2;
	//	}
	//	return coordinate[0] / (m_shape[0] / (xDivisions + 1)) + (coordinate[1] / (m_shape[1] / (yDivisions + 1))) + xDivisions * (coordinate[0] / (m_shape[0] / (xDivisions + 1)));
	//}
	//return 0;

	//if (totalIndex(coordinate) == m_totalVolume) {
	//	return mpiWrapper::nProcs();
	//}

	int out = 0;
	int temp;
	for (int i = 0; i < m_Ndims; i++) {
		temp = coordinate[i] / (m_shape[i] / (m_cuts[i] + 1));
		for (int j = 0; j < i; j++) {
			temp *=(m_cuts[j] + 1);
		}
		out += temp;
	}
	//if (mpiWrapper::id() == 0) {
	//	std::cout << coordinate[0] << " " << coordinate[1] << " " << coordinate[2] << " " << coordinate[3] << "		:		" << out << std::endl;
	//}
	
	return out;
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
	if (type == "torus") {
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
	if (type == "torus_OpenT") {
		if (totalIndex(coordinate) == m_totalVolume) {
			fwd_coor[m_Ndims - 1] = m_shape[m_Ndims - 1];

			for (int i = 0; i < m_Ndims - 1; i++) {
				fwd_coor[i] = m_shape[i] - 1;
			}
			back_coor[m_Ndims - 1] = m_shape[m_Ndims - 1];
			for (int i = 0; i < m_Ndims - 1; i++) {
				back_coor[i] = m_shape[i] - 1;
			}
		}
		if (0 == direction) {
			if (coordinate[0]+1 == m_shape[0]) {
				
				fwd_coor[m_Ndims-1] = m_shape[m_Ndims-1];
				
				for (int i = 0; i < m_Ndims-1; i++) {
					fwd_coor[i] = m_shape[i]-1;
				}

				//fwd_coor[0] = m_shape[0];
				//for (int i = 1; i < m_Ndims; i++) {
				//	fwd_coor[i] = m_shape[i] - 1;
				//}

			}
			else {
				fwd_coor[0] = coordinate[0] + 1;
			}
			if (coordinate[0] == 0) {
				back_coor[m_Ndims-1] = m_shape[m_Ndims-1];
				for (int i = 0; i < m_Ndims-1; i++) {
					back_coor[i] = m_shape[i] - 1;
				}

				//back_coor[0] = m_shape[0];
				//for (int i = 1; i < m_Ndims; i++) {
				//	back_coor[i] = m_shape[i] - 1;
				//}


				//std::cout << "coor :" << coordinate[0] << " " << coordinate[1] << " " << coordinate[2] << " " << coordinate[3] << std::endl;
				//std::cout << "back_coor :" << back_coor[0] << " " << back_coor[1] << " " << back_coor[2] << " " << back_coor[3] << std::endl;
				//std::cout << "total Idx of back_coor :	" << totalIndex(back_coor) << std::endl;


			}
			else {
				back_coor[0] = coordinate[0] - 1;
			}
		}
		else {
			fwd_coor[0] = coordinate[0];
			back_coor[0] = coordinate[0];
			for (int i = 1; i < m_Ndims; i++) {
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

	}

}

bool Lattice::is_SharedMemory(int* coordinate, int* coorFwd, int* coorBack){
	int id_fwd;
	int id_back;
	bool out = false;
	if (totalIndex(coordinate) == m_totalVolume) {
		return false;
	}

	for (int i = 0; i < m_Ndims; i++) {
		NearestNeighbour(coordinate, i, coorFwd, coorBack);
		id_fwd = Coordinate_ProcID(coorFwd);
		id_back = Coordinate_ProcID(coorBack);
		if (totalIndex(coorFwd) == m_totalVolume || totalIndex(coorBack) == m_totalVolume) {
			out = false;
		}
		if (mpiWrapper::id() == id_fwd || mpiWrapper::id() == id_back) {
			out = true;
		}
				for (int j = 0; j < m_Ndims; j++) {
					if (j != i || m_cloversRequired) {
					NearestNeighbour(coorFwd, j, m_coorFwdFwd, m_coorFwdBack);
					NearestNeighbour(coorBack, j, m_coorBackFwd, m_coorBackBack);
					if (totalIndex(m_coorFwdFwd) == m_totalVolume ||
						totalIndex(m_coorBackBack) == m_totalVolume ||
						totalIndex(m_coorFwdBack) == m_totalVolume ||
						totalIndex(m_coorBackFwd) == m_totalVolume) {
						out = false;
					}
					if (Coordinate_ProcID(m_coorFwdFwd) == mpiWrapper::id() ||
						Coordinate_ProcID(m_coorBackBack) == mpiWrapper::id() ||
						Coordinate_ProcID(m_coorFwdBack) == mpiWrapper::id() ||
						Coordinate_ProcID(m_coorBackFwd) == mpiWrapper::id()) {
						out = true;
					}
				}
					if (m_cloversRequired) {
					for (int k = 0; k < m_Ndims; k++) {
						NearestNeighbour(m_coorBackFwd, k, m_coorBackFwdFwd, m_coorBackFwdBack);
						NearestNeighbour(m_coorFwdFwd, k, m_coorFwdFwdFwd, m_coorFwdFwdBack);
						NearestNeighbour(m_coorBackBack, k, m_coorBackBackFwd, m_coorBackBackBack);
						if (totalIndex(m_coorFwdFwdFwd)		== m_totalVolume ||
							totalIndex(m_coorFwdFwdBack)	== m_totalVolume ||
							totalIndex(m_coorBackBackFwd)	== m_totalVolume ||
							totalIndex(m_coorBackBackBack)	== m_totalVolume ||
							totalIndex(m_coorBackFwdFwd)	== m_totalVolume ||
							totalIndex(m_coorBackFwdBack)	== m_totalVolume) {
							out = false;
						}
						if (Coordinate_ProcID(m_coorFwdFwdFwd)		== mpiWrapper::id() ||
							Coordinate_ProcID(m_coorFwdFwdBack)		== mpiWrapper::id() ||
							Coordinate_ProcID(m_coorBackBackFwd)	== mpiWrapper::id() ||
							Coordinate_ProcID(m_coorBackBackBack)	== mpiWrapper::id() ||
							Coordinate_ProcID(m_coorBackFwdFwd)		== mpiWrapper::id() ||
							Coordinate_ProcID(m_coorBackFwdBack)	== mpiWrapper::id() ) {
							out = true;
						}
							for (int l = 0; l < m_Ndims; l++) {
								if (l != k) {
									NearestNeighbour(m_coorFwdFwdFwd, l, m_coorFwdFwdFwdFwd, m_coorFwdFwdFwdBack);
									NearestNeighbour(m_coorBackBackBack, l, m_coorBackBackBackFwd, m_coorBackBackBackBack);
									NearestNeighbour(m_coorBackFwdBack, l, m_coorBackFwdBackFwd, m_coorBackFwdBackBack);
									if (totalIndex(m_coorFwdFwdFwdFwd)		== m_totalVolume ||
										totalIndex(m_coorBackBackBackBack)	== m_totalVolume ||
										totalIndex(m_coorFwdFwdFwdBack)		== m_totalVolume ||
										totalIndex(m_coorBackFwdBackFwd)	== m_totalVolume ||
										totalIndex(m_coorBackBackBackFwd)	== m_totalVolume) {
										out = false;
									}
									if (Coordinate_ProcID(m_coorFwdFwdFwdFwd)		== mpiWrapper::id() ||
									/*	Coordinate_ProcID(m_coorBackFwdBackBack)	== mpiWrapper::id() ||*/
										Coordinate_ProcID(m_coorBackBackBackBack)	== mpiWrapper::id() ||
										Coordinate_ProcID(m_coorFwdFwdFwdBack)		== mpiWrapper::id() ||
										Coordinate_ProcID(m_coorBackFwdBackFwd)		== mpiWrapper::id() ||
										Coordinate_ProcID(m_coorBackBackBackFwd)	== mpiWrapper::id()) {
										out = true;
									}
								}
							}
					}
			}
		}
	}
	return out;
}


int* Lattice::getShape() const{
	return m_shape;
}

int Lattice::getNdims() const{
	return m_Ndims;
}

//int Lattice::getthisProc_Volume() const{
//	return m_thisProc_Volume;
//}

std::string Lattice::getType() const{
	return type;
}

void Lattice::Print_Partitioning() {
	if (mpiWrapper::id() == 0) {
		int* coordinate = new int[4];
		coordinate[2] = 0;
		coordinate[3] = 0;
		for (int i = 0; i < m_shape[0]; i++) {
			for (int j = 0; j < m_shape[1]; j++) {
				coordinate[0] = i;
				coordinate[1] = j;
				std::cout << Coordinate_ProcID(coordinate) << " ";
			}
			std::cout << "\n";
		}
		delete[] coordinate;
		std::cout.flush();
	}
	MPI_Barrier(mpiWrapper::comm());
}