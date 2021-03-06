#pragma once 
#include "../mpiWrapper.h"
#include "../SimulationParameters.h"
#include <string>


class Lattice {
public:
	Lattice(SimulationParameters &params, bool cloversRequired = true);

	void partition_lattice();
	//Tells each process which sites (in its own internal index) are shared with which neighbouring processes
	//and stores it in the array  m_Buffer_receive by proc_id and a index from zero to Nr. of sites shared with that process
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
	void Print_Partitioning();
	int* getShape() const;
	int getNdims() const;
	//int getthisProc_Volume() const;
	inline int getthisProc_Volume() const{
	return m_thisProc_Volume;
	}
	inline int getthisProc_RespVolume() const {
		return m_responsible_Volume;
	}
	std::string getType() const;
	

	int m_totalVolume;//The total number of lattice points for the entire lattice
	int m_thisProc_Volume;//The number of lattice points stored by this process, including shared points
	int m_responsible_Volume;//The number of lattice points stored by this process, excluding shared points
	int* m_thisProc_TotalIndex;
	int* m_InternalToTotal_idx;//gives the total index corresponding to that internal index
	int* m_TotalToInternal_idx;//gives the internal index corresponding to that total index (assuming one exists)
	int* m_Internal_IdxToProcId;//gives the processor identity of the lattice point with that internal index
	int* m_Internal_IdxToParity;//gives the parity of the lattice point with that internal index
	
	int** m_cntr;//The center coordinate by [Internal index][space-time direction]
	
	//returns the internal index of the 
	int** m_fwd;//the forward neighbour by [Internal index][space-time direction]
	int** m_back;//the backward neighbour by [Internal index][space-time direction]
	// moving forward or backwards in the space-time direction given

	int** m_InternalIdx_start;//[procID][parity(=0/1)]->the starting internal index for boundary sites governed by procID with that parity
	int** m_InternalIdx_stop;//[procID][parity(=0/1)]->the last internal index for boundary sites governed by procID with that parity


	int* m_BufferSize;//The size of the buffer that should be allocated for sending/receiving for each process, including itself and processes that don't share sites ID i.e. size=0
	/// <summary>
	/// indexed by [proc id][i] where i runs from 0 to m_BufferSize[proc id]:
	/// Tells THIS process, which sites <proc id> needs to receive information about in THIS process' internal index
	/// 
	/// USE:
	/// for i=0 -> i< m_buffersize[sending_to_id] -> i++ :
	///		send field value at site m_Buffer_receive[sending_to_id][i]
	/// </summary>
	int** m_Buffer_receive;


protected:
	bool m_cloversRequired;
	//bool Trivial_Flag;
	bool TPancake_Flag;
	//bool STFF_Flag;
	void CheckPartitioning();
	int* m_cuts;
	int* m_shape;//The size of the lattice for each dimension
	
	int m_Ndims;//The number of lattice dimensions

	int* m_coordinate;//stores the current working coordinate 
	int* m_coorFwd;//stores a coordinate forwardly neighbouring the current working coordinate 
	int* m_coorBack;//stores a coordinate backwardly neighbouring the current working coordinate 
	int* m_coorFwdFwd;//stores a coordinate forwardly neighbouring the current working coordinate 
	int* m_coorBackBack;//stores a coordinate backwardly neighbouring the current working coordinate 
	int* m_coorFwdBack;//stores a coordinate forwardly neighbouring the current working coordinate 
	int* m_coorBackFwd;//stores a coordinate backwardly neighbouring the current working coordinate 
	int* m_coorFwdFwdFwd;
	int* m_coorBackBackBack;
	int* m_coorFwdFwdBack;
	int* m_coorBackBackFwd;

	int* m_coorBackFwdFwd;
	int* m_coorBackFwdBack;
	int* m_coorBackFwdBackFwd;
	int* m_coorFwdFwdFwdFwd;//NOT IN USE
	int* m_coorBackFwdBackBack;//NOT IN USE
	int* m_coorBackBackBackBack;
	int* m_coorFwdFwdFwdBack;
	int* m_coorBackBackBackFwd;

	std::string type;




};