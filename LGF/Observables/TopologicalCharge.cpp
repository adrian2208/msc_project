#include "TopologicalCharge.h"


TopologicalCharge::TopologicalCharge(SU3_field& U) {
	m_AutoCorrTime = calculate_AutoCorrTime();
	m_U = &U;
}

void TopologicalCharge::calculate() {
	static const double PI = 4.0 * atan(1.0);
	static const double PI2 = PI * PI;
	double localSum = 0.0;
	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
		localSum += ((*m_U).clover_avg(i, 0, 1) * (*m_U).clover_avg(i, 2, 3)).ReTr();
		localSum -= ((*m_U).clover_avg(i, 0, 2) * (*m_U).clover_avg(i, 1, 3)).ReTr();
		localSum += ((*m_U).clover_avg(i, 0, 3) * (*m_U).clover_avg(i, 1, 2)).ReTr();
	}
	double totalSum;
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	totalSum *= 8.0 / (32.0 * PI2);
	m_resultVector.push_back(totalSum);
	if (mpiWrapper::id() == 0) {
		std::cout << totalSum << "\n";
	}
	//double avg_Plaquette = (*m_U).Avg_Plaquette();
	//std::cout << "exp_val_Plaquette: " << avg_Plaquette << "\n\n";
}

int TopologicalCharge::calculate_AutoCorrTime(){
	return 0;
}
