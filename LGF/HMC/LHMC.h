#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Field/Double_field.h"
#include "../Action/Action.h"
#include "../Action/Gauge/Wilson.h"
#include "../Math/SU3_gen.h"
class LHMC {
public:
	LHMC(SU3_field& U, Wilson& GaugeAction, double epsilon,int NrLeaps = 0);
	~LHMC();
	//void Initialize_P();
	double calculate_KineticEnergy(su3_mat& P);
	void leapfrog(SU3_field& U, su3_mat& P, su3_mat& F, Wilson& action, int parity, int i, int mu);
	void update(int parity);
	void sweep();
	double acceptanceRate();
private:
	SU3_field* m_U;
	//RNG_field* m_rand; I would use this, but it takes up too much RAM
	Random* m_rand;
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