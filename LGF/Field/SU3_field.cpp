#include "SU3_field.h"

SU3_field::SU3_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF){}

void SU3_field::saveSU3ToFile(const std::filesystem::path& identifier){
	std::filesystem::path fieldType("/su3/");
	saveToFile(fieldType, identifier);
}

void SU3_field::loadSU3FromFile(const std::filesystem::path& identifier){
	std::filesystem::path fieldType("/su3/");
	loadFromFile(fieldType, identifier);
}
