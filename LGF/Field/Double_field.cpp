#include "Double_field.h"

Double_field::Double_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF) {
	for (int i = (*this).Responsible_Start(); i < (*this).Responsible_Stop(); i++) {
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			(*this)(i, mu) = 0.0;
		}
	}
}

void Double_field::operator=(const Double_field& field) {
	for (int i = 0; i < field.m_NrExtDOF * field.m_lattice->getthisProc_Volume(); i++) {
		FieldArray[i] = field.FieldArray[i];
	}
}
