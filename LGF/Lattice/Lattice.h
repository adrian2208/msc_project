#pragma once 
#include "../mpiWrapper.h"



class Lattice {
public:
	Lattice(int Ndims, int shape[]);

	void partition_lattice();
	//Assigns a single index based on the coordinates according to scheme outlined in function definition
	int totalIndex(int *coordinate);
	
	//This virtual method allocates processes to coordinates according to the first coordinate
	virtual int Coordinate_ProcID(int* coordinate);
	//This virtual method encodes periodic boundary conditions for all directions on the lattice
	virtual void NearestNeighbour(int* coordinate, int direction, int* fwd_coor, int* back_coor);



	void print_indices();

	void print();
	
	int* m_shape;//The size of the lattice for each dimension
	int* m_coordinate;//stores the current working coordinate 
	int m_Ndims;//The number of lattice dimensions
	int m_totalVolume;//The total number of lattice points for the entire lattice
	int m_thisProc_Volume;//The number of lattice points per subprocess
	int* m_thisProc_TotalIndex;

	int** m_fwd;//the forward neighbour by [Process_lattice_index][space-time_direction]
	int** m_back;//the backward neighbour by [Process_lattice_index][space-time_direction]

};