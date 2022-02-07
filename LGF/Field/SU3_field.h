#pragma once
#include "Field.h"
#include "../Math/SU3_mat.h"
#include "../Math/SU3_gen.h"
#include "../Math/Random.h"
#include <vector>
class SU3_field : public Field<su3_mat> {
public:
	SU3_field(Lattice& lattice, int NrExtDOF);//NrExtDOF corresponds to e.g. mu={0,1,2,3}
	void saveSU3ToFile(double beta, const std::string& identifier = "", const std::string& dataFolder = "C:/Users/adria/Documents/msc_project/data/");
	void loadSU3FromFile(double beta, const std::string& identifier = "", const std::string& dataFolder = "C:/Users/adria/Documents/msc_project/data/");
	void transfer_FieldValues();
	void InitializeHotStart();
	void InitializeColdStart();

	/// <summary>
	/// IMPORTANT: This function will not update the boundary field values before calculating!
	/// Make sure that the boundary field values are up to date before calling this function.
	/// </summary>
	double total_PlaquetteSum();
	double Avg_Plaquette();
	su3_mat staple(int internal_index, int mu);
	su3_mat clover_avg(int internal_index, int mu, int nu);
	//su3_mat clover_avg1(int internal_index, int mu, int nu);
	inline su3_mat plaquette(int internal_index, int mu, int nu);
	
	void operator= (const SU3_field& field);


};