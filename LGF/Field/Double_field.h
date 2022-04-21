#pragma once
#include "Field.h"
class Double_field : public Field<double> {
public:
	Double_field(Lattice& lattice, int NrExtDOF);
	void saveDoubleToFile(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier);
	void operator= (const Double_field& field);

};