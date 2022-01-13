#pragma once
#include "../Action.h"
#include "../../Field/SU3_field.h"
#include "../../Math/SU3_mat.h"
class Wilson : public Action {
public:
	Wilson(double beta);
	/// <summary>
	/// If the program is getting bad values:
	/// reasses the coefficients here as well as those included by the total_PlaquetteSum() function
	/// </summary>
	/// <returns>the Action</returns>
	double calculate_Action(SU3_field& field);
	void calculate_Force(SU3_field& field, SU3_field& F);
private:
	double m_beta;
};