#include "RNG_field.h"

RNG_field::RNG_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF) {
	std::random_device seed;
	for (int i = (*this).Responsible_Start(); i < (*this).Responsible_Stop(); i++){
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			(*this)(i, mu) = Random(seed() + mpiWrapper::id()+i * m_NrExtDOF + mu);
		}
	}
}
