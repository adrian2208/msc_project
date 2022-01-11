#pragma once
#include "../Field/SU3_field.h"
#include "../Action/Action.h"
class HMC {
public:
	HMC(SU3_field& U,Action& GaugeAction, double LFStepSize);
	void Initialize_P;

private:
	SU3_field* P;
	SU3_field* F;

	int m_NrAccepted;
	int m_NrstepsTaken;
	double m_LFStepSize
};