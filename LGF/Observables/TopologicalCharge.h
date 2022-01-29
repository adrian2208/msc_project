#pragma once
#include "../Field/SU3_field.h"
class TopologicalCharge {

public:
	TopologicalCharge(SU3_field& U);
	void calculate();
	int calculate_AutoCorrTime();

private:
	SU3_field* m_U;
	int m_AutoCorrTime;
	std::vector<double> m_resultVector;
};