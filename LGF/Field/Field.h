#pragma once
#include "../Lattice/Lattice.h"
#include "../Vertex.h"
template<class T>
class Field {
public:
	Field(Lattice& lattice, int NrExtDOF) {
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
	int m_NrExtDOF;// e.g. 4: mu=0,1,2,3
	int size_T;//the number of real variables required to store T
	int size_localField;


};