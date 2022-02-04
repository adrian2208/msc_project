#pragma once
#include "../Field/SU3_field.h"
class EnergyDensity {

public:
	EnergyDensity(SU3_field& U);
	void calculate(double flowTime);
	int calculate_AutoCorrTime();
	void saveEnergyDensityToFile(double beta, const std::filesystem::path& identifier = "");
private:
	SU3_field* m_U;
	int m_AutoCorrTime;
	std::vector<double> m_resultVector;
};