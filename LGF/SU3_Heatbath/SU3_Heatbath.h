#pragma once
#include "../Field/SU3_field.h"
#include "../Field/RNG_field.h"
#include "../Field/Double_field.h"
#include "../Action/Action.h"
#include "../Action/Gauge/Wilson.h"
#include "../Math/SU3_gen.h"
class SU3_Heatbath {
public:
	SU3_Heatbath(SU3_field& U);

	void update(SU3_field& U, int n_iter, double beta);
	void heatbath_SU2(Random& random, double beta_eff, C_double* a);

private:
	RNG_field* m_rand;
};
