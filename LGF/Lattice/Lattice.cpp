#include "Lattice.h"

Lattice::Lattice(int Ndims, int shape[]) {
	
	//Storing local variables
	m_Ndims = Ndims;
	m_shape = new int[Ndims];
	m_coordinate = new int[Ndims];
	//Calculates the total volume of the lattice & writing to member variables
	m_totalVolume = 1;
	m_thisProc_Volume = 0;
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
		if (mpiWrapper::id() == id) {							//  |
			m_thisProc_TotalIndex[m_thisProc_Volume] = totalIdx;//  |-- identifies the process working on this
			m_thisProc_Volume++;								//  |	coordinate, incrementing the proc volume.
		}														//__|
		//returns new m_coordinate by reference eventually traversing the lattice															
		traverse_lattice(m_coordinate, m_Ndims, m_shape); 												
		i++;
	}

	//Initializing arrays for storing nearest neighbour information
	m_fwd = new int* [m_thisProc_Volume];
	m_back = new int* [m_thisProc_Volume];
	for (i = 0; i < m_thisProc_Volume; i++) {
		m_fwd[i] = new int[m_Ndims];
		m_back[i] = new int[m_Ndims];
	}
	
	for (int proc_id = 0; proc_id < mpiWrapper::nProcs(); proc_id++) {
		for (int par = 0; par <= 1; par++) {

		}
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
