#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Field/Double_field.h"
#include "../Action/Action.h"
#include "../Action/Gauge/Wilson.h"
#include "../Math/SU3_gen.h"
class SU3_Heatbath {
public:
	SU3_Heatbath(SU3_field& U,double beta,int OR_per_HB);

	void update(int n_iter);
	void Cabibbo_Marinari(SU3_field& U);
	void OverRelaxation(SU3_field& U);
	void heatbath_SU2(C_double* TSR, int internIdx);

	std::string getupdateMethod() const;

private:
	std::string m_updateMethod;
	int m_OR_per_HB;
	double m_beta;
	SU3_field* m_U;
	RNG_field* m_rand;
};
