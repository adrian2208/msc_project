#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Field/Double_field.h"
#include "../Action/Action.h"
#include "../Action/Gauge/Wilson.h"
#include "../Math/SU3_gen.h"
class LHMC {
public:
	LHMC(SU3_field& U, Wilson& GaugeAction, double epsilon);
	//void Initialize_P();
	double calculate_KineticEnergy(SU3_field& P, int i, int mu);
	void leapfrog(SU3_field& U, SU3_field& P, SU3_field& F, Wilson& action, int parity, int i, int mu);
	void update(int parity);
	void sweep();
	double acceptanceRate();
private:
	SU3_field* m_P;
	SU3_field* m_F;
	SU3_field* m_U;
	SU3_field* m_F_init;
	SU3_field* m_U_init;
	RNG_field* m_rand;
	Double_field* m_K_init;
	Double_field* m_localEpsilon;
	//Double_field* m_S_init;
	Wilson* m_GaugeAction;

	int m_NrAccepted;
	int m_NrstepsTaken;
	int m_NrleapfrogSteps;
	double m_epsilon;



	double avg_Plaquette;
#ifdef _DEBUG

	double exp_val_expdeltaH;
#endif // _DEBUG

};