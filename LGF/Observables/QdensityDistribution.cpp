#include "QdensityDistribution.h"
#include <fstream>

QdensityDistribution::QdensityDistribution(SU3_field& U,double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {
	m_updateMethod = updateMethod;
	m_identifier = identifier;
	m_dataFolder = dataFolder;
	m_U = &U;
	m_beta = beta;
	m_q = new Double_field(U.getLatticePtr(), 1);
}
QdensityDistribution::~QdensityDistribution() {
	delete m_q;
}
void QdensityDistribution::calculate(double flowTime) {
	static const double PI = 4.0 * atan(1.0);
	static const double PI2 = PI * PI;
	double localSum;
	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) { //This loop, alongside the multiple of eight constitutes the 4D levi cevita, thanks to the cyclicity of the trace 
		localSum = 0.0;
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 2)) * ((*m_U).Improved_fieldStrengthTensor(i, 0, 1))).ReTr();
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 1)) * ((*m_U).Improved_fieldStrengthTensor(i, 2, 0))).ReTr();
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 0)) * ((*m_U).Improved_fieldStrengthTensor(i, 1, 2))).ReTr();
		(*m_q)(i, 0) = (1.0/4.0*PI2)*localSum;
	}

	(*m_q).saveDoubleToFile(m_beta, m_updateMethod, m_dataFolder, m_identifier);
}