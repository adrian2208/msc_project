#pragma once
#include "Lattice/Lattice.h"
//vertex functionality:
// operator +, - etc.

class Vertex {
public:
	Vertex(const Lattice& lattice);
	int index;
private:
	Lattice* m_lattice_ptr;
};