#pragma once
#include "Field.h"
class Double_field : public Field<double> {
public:
	Double_field(Lattice& lattice, int NrExtDOF);

	void operator= (const Double_field& field);

};