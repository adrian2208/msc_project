#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Action/Action.h"
#include "../Math/SU3_gen.h"
class HMC {
public:
	HMC(SU3_field& U,Action& GaugeAction, double epsilon);
	//void Initialize_P();
	double calculate_KineticEnergy(SU3_field& P);
	void leapfrog(SU3_field& U, SU3_field& P,SU3_field& F, Action& action);
	void update();
	double acceptanceRate();
private:
	SU3_field* m_P;
	SU3_field* m_F;
	SU3_field* m_U;
	SU3_field* m_F_init;
	SU3_field* m_U_init;
	RNG_field* m_rand;
	Action* m_GaugeAction;

	double S_init;

	int m_NrAccepted;
	int m_NrstepsTaken;
	int m_NrleapfrogSteps;
	double m_epsilon;



	double exp_val_Plaquette;
#ifdef _DEBUG
	
	double exp_val_expdeltaH;
#endif // _DEBUG

};