#include "RNG_field.h"

RNG_field::RNG_field(Lattice& lattice) : Field(lattice, 1) {
	std::random_device seed;
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		(*this)(i, 0) = Random(seed()+i);
	}
}
