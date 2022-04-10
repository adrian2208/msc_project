#pragma once
#include "../Field/SU3_field.h"
#include "../Field/Double_field.h"
class QdensityDistribution {

public:
	QdensityDistribution(SU3_field& U, double beta, const std::string& updateMethod, const std::string& identifier = "", const std::string& dataFolder = "C:/Users/adria/Documents/msc_project/data/");
	~QdensityDistribution();
	void calculate(double flowTime);
private:
	SU3_field* m_U;
	Double_field* m_q;
	double m_beta;
	std::string m_updateMethod;
	std::string m_identifier;
	std::string m_dataFolder;
};