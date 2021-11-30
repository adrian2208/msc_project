#include "Lattice.h"

Lattice::Lattice(int Ndims, int shape[]) {
	
	//Storing local variables
	m_Ndims = Ndims;
	m_shape = new int[Ndims];
	m_coordinate = new int[Ndims];
	//Calculates the total volume of the lattice & writing to member variables
	m_totalVolume = 1;
	for (int i = 0; i < Ndims; i++) {
		m_shape[i] = shape[i];
		m_totalVolume *= shape[i];
		m_coordinate[i] = 0;
	}
}

void Lattice::partition_lattice(int NrPartitions){
	long int* thisProc_TotalIndex = new long int[m_totalVolume];
	int totalIdx;
	m_thisProc_Volume = 0;

	for (int i = 0; i < m_Ndims; i++) {
		for (int j = 0; j < m_shape[i]; j++) {
			m_coordinate[i] = j;
			totalIdx = totalIndex(m_coordinate);
			
			if (mpiWrapper::id() == Coordinate_ProcID(m_coordinate)) {
				thisProc_TotalIndex[m_thisProc_Volume] = totalIdx;
				m_thisProc_Volume++;
			}
		}
		m_coordinate[i] = 0;
	}
	m_fwd = new int* [m_thisProc_Volume];
	m_back = new int* [m_thisProc_Volume];
	for (int i = 0; i < m_thisProc_Volume; i++) {
		m_fwd[i] = new int[m_Ndims];
		m_back[i] = new int[m_Ndims];
	}

}

int Lattice::Coordinate_ProcID(int* coordinate){
	
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

void Lattice::NearestNeighbour(){
}


