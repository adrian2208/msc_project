#include "Vertex.h"

Vertex::Vertex(const Lattice& lattice){
	m_lattice_ptr = (Lattice*) &lattice;
	index = (*m_lattice_ptr).start[ME][0];
}
