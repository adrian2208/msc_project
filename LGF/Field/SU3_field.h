#include "Field.h"
#include "../Math/SU3_mat.h"
#pragma once
class SU3_field : public Field<su3_mat> {
public:
	SU3_field(Lattice& lattice, int NrExtDOF);//NrExtDOF corresponds to e.g. mu={0,1,2,3}
	void saveSU3ToFile(const std::filesystem::path& identifier = "");
	void loadSU3FromFile(const std::filesystem::path& identifier = "");
	void CommunicateFieldValues();

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!TRY WITHOUT MPI_ALLREDUCE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double total_PlaquetteSum();
	double Avg_Plaquette();
	su3_mat staple(int internal_index, int mu);
	inline su3_mat plaquette(int internal_index, int mu, int nu);
	
	void operator= (const SU3_field& field);
	//void  copyFieldVals(const SU3_field& field);
};