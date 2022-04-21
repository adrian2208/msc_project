#pragma once
#include "../Field/SU3_field.h"
class TopologicalCharge {

public:
	TopologicalCharge(SU3_field& U);
	void calculate(double flowTime);
	void saveTopologicalChargeToFile(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier="");
private:
	SU3_field* m_U;
	
	std::vector<double> m_resultVector;
	std::vector<double> m_FlowMeasurementTimeVector;
};