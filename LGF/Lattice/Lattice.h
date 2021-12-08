#pragma once 
#include "../mpiWrapper.h"



class Lattice {
public:
	Lattice(int Ndims, int shape[]);

	void partition_lattice();
	void DistributeBuffers();
	//Assigns a single index based on the coordinates according to scheme outlined in function definition
	int totalIndex(int *coordinate);
	//returns the total coordinate by reference from the total index
	void IndexToCoordinate(int index, int* coordinate);
	
	//This virtual method allocates processes to coordinates according to the first coordinate
	virtual int Coordinate_ProcID(int* coordinate);
	//This virtual method encodes periodic boundary conditions for all directions on the lattice
	virtual void NearestNeighbour(int* coordinate, int direction, int* fwd_coor, int* back_coor);
	//This virtual method sets the size of the overlap stored by neighbouring processes
	// defaulting to 1 vertex overlap
	virtual bool is_SharedMemory(int* coordinate, int* coorFwd, int* coorBack);



	void print_indices();

	void print();
	
	int* m_shape;//The size of the lattice for each dimension
	int m_Ndims;//The number of lattice dimensions
	int m_totalVolume;//The total number of lattice points for the entire lattice
	int m_thisProc_Volume;//The number of lattice points stored by this process, including shared points
	int m_responsible_Volume;//The number of lattice points stored by this process, excluding shared points
	int* m_thisProc_TotalIndex;
	int* m_InternalToTotal_idx;//gives the total index corresponding to that internal index
	int* m_TotalToInternal_idx;//gives the internal index corresponding to that total index (assuming one exists)
	int* m_Internal_IdxToProcId;//gives the processor identity of the lattice point with that internal index
	int* m_Internal_IdxToParity;//gives the parity of the lattice point with that internal index
	
	int** m_cntr;//The center coordinate by [Internal index][space-time direction]
	int** m_fwd;//the forward neighbour by [Internal index][space-time direction]
	int** m_back;//the backward neighbour by [Internal index][space-time direction]

	int** m_InternalIdx_start;
	int** m_InternalIdx_stop;
	int** m_BufferSize;
	int** m_Buffer_receive;


protected:
	int* m_coordinate;//stores the current working coordinate 
	int* m_coorFwd;//stores a coordinate forwardly neighbouring the current working coordinate 
	int* m_coorBack;//stores a coordinate backwardly neighbouring the current working coordinate 

};