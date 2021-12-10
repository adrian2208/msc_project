#pragma once
#include "../Lattice/Lattice.h"
#include "../Vertex.h"
template<class T>
class Field {
public:
	Field(Lattice& lattice, int NrExtDOF) {
		m_lattice = &lattice;
		m_NrExtDOF = NrExtDOF;
		allocate_field(lattice, NrExtDOF);
	}
	void allocate_field(Lattice& lattice, int NrExtDOF) {
		size_localField =  NrExtDOF* lattice.m_thisProc_Volume;
		container = new T[size_localField];
	}
	
	// overloads () operator on Field instance to return T situated on Vertex x and ext. DOF mu by reference
	inline T& operator() (Vertex x, int mu) {
		return container[x.index * m_NrExtDOF + mu];
	}

protected:
	T* container;
	int m_NrExtDOF;//the number of external DOF e.g. 4: mu=0,1,2,3
	//the number of internal DOF e.g. 1 for an integer field, 2 for a complex field or 
	//2*9=18 for an SU(3) field in the redundant representation (otherwise 3^2-1=8, but this is not as efficient)
	int m_NrIntDOF;
	int m_sizeOfT;//The number of bytes needed to store the type T
private:
	Lattice* m_lattice;

};