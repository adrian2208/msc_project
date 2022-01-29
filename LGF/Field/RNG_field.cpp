#include "RNG_field.h"

RNG_field::RNG_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF) {
	std::random_device seed;
	for (int i = 0; i < lattice.m_thisProc_Volume; i++) {
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			(*this)(i, mu) = Random(seed() + i * m_NrExtDOF + mu);
		}
	}
}
