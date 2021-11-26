#pragma once
#include "mpi.h"

class Lattice {
public:
	Lattice(int Ndims, int shape[]);

	virtual void partition_lattice(int NrPartitions);
	int Coordinate_ProcID(int* coordinate);//returns ID of processor responsible for this coordinate
	int totalIndex(int *coordinate);//Assigns a single index based on the coordinates according to scheme outlined in function definition
	
	virtual void NearestNeighbour() = 0;//Encodes boundary conditions
	
	int* m_shape;
	int* m_coordinate;
	int m_Ndims;
	int m_totalVolume;//The total number of lattice points for the entire lattice
	int m_thisProc_Volume;//The number of lattice points per subprocess

};