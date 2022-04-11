#pragma once
#include "../Field/SU3_field.h"
class ImprovedTopologicalCharge {

public:
	ImprovedTopologicalCharge(SU3_field& U);
	void calculate(double flowTime);
	int calculate_AutoCorrTime();
	void saveTopologicalChargeToFile(double beta, const std::string& updateMethod, const std::string& identifier = "", const std::string& dataFolder = "C:/Users/adria/Documents/msc_project/data/");
private:
	SU3_field* m_U;
	int m_AutoCorrTime;

	std::vector<double> m_resultVector;
};