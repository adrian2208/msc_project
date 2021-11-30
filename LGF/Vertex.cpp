#include "Vertex.h"

Vertex::Vertex(const Lattice& lattice){
	m_lattice_ptr = (Lattice*) &lattice;
	index = 0;//(*m_lattice_ptr).start[ME][0];
}
