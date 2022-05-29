#include "FastGF.h"

/// 
/// IN DEVELOPMENT!!!
/// DO NOT USE!
/// 


FastGF::FastGF(Action& GaugeAction, SU3_field& U, double epsilon) {
	m_updateMethod = "GF";
	m_flowTime = 0.0;
	m_NrSteps = 0;
	m_GaugeAction = &GaugeAction;
	m_epsilon = epsilon;
	m_U = &U;
	m_Z0 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_Z1 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_Z2 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	(*m_U).transfer_FieldValues();
}

void FastGF::flow() {

	su3_mat temp;
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int parity = 0; parity < 2; parity++) {
			(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z0), parity, mu);
			for (int i = (*m_U).Responsible_Start(parity); i < (*m_U).Responsible_Stop(parity); i++) {
				temp = 0.25 * m_epsilon * (*m_Z0)(i, mu);
				(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);

			}
		}
	}

	(*m_U).transfer_FieldValues();
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int parity = 0; parity < 2; parity++) {
			(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z1), parity, mu);
			for (int i = (*m_U).Responsible_Start(parity); i < (*m_U).Responsible_Stop(parity); i++) {
				temp = m_epsilon * (8.0 / 9.0 * (*m_Z1)(i, mu) - 17.0 / 36.0 * (*m_Z0)(i, mu));
				(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);

			}
		}
	}

	
	(*m_U).transfer_FieldValues();
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int parity = 0; parity < 2; parity++) {
			(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z2), parity, mu);
			for (int i = (*m_U).Responsible_Start(parity); i < (*m_U).Responsible_Stop(parity); i++) {
				temp = m_epsilon * (3.0 / 4.0 * (*m_Z2)(i, mu) - 8.0 / 9.0 * (*m_Z1)(i, mu) + 17.0 / 36.0 * (*m_Z0)(i, mu));
				(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);

			}
		}
	}

	m_flowTime += m_epsilon;
	m_NrSteps++;
	(*m_U).transfer_FieldValues();

	MakeMeasurements();
}
double FastGF::GetFlowTime() const {
	return m_flowTime;
}

void FastGF::Include_TopCharge(TopologicalCharge& TopCharge) {
	(*m_U).transfer_FieldValues();
	m_topologicalCharge = &TopCharge;
}
void FastGF::Include_EnergyDensity(EnergyDensity& Edensity) {
	(*m_U).transfer_FieldValues();
	m_EnergyDensity = &Edensity;
}

void FastGF::MakeMeasurements() {
	if (m_topologicalCharge != nullptr) {
		m_topologicalCharge->calculate(m_flowTime);
	}
	if (m_EnergyDensity != nullptr) {
		m_EnergyDensity->calculate(m_flowTime);
	}
}
std::string FastGF::getupdateMethod() const {
	return m_updateMethod;
}