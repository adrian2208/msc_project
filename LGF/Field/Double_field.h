#pragma once
#include "Field.h"
class Double_field : public Field<double> {
public:
	Double_field(Lattice& lattice, int NrExtDOF);
	void saveDoubleToFile(double beta, const std::string& updateMethod, const std::string& identifier, const std::string& dataFolder);
	void operator= (const Double_field& field);

};