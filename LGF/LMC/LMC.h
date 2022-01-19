#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Action/Action.h"
#include "../Math/SU3_gen.h"
class LMC {
public:
	LMC(SU3_field& U, double epsilon);
	void update();
	double acceptanceRate();
private:
	
	SU3_field* m_U;

	int m_NrAccepted;
	int m_NrstepsTaken;
	double m_epsilon;


	double exp_val_Plaquette;
	double exp_val_expdeltaH;


};