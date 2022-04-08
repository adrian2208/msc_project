#include "GradientFlow.h"

GradientFlow::GradientFlow(Action& GaugeAction, SU3_field& U, double epsilon, int FlowStepsPerMeasurement){
	m_updateMethod = "GF";
	m_FlowStepsPerMeasurement = FlowStepsPerMeasurement;
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
GradientFlow::~GradientFlow() {
	delete m_Z0;
	delete m_Z1;
	delete m_Z2;
}

void GradientFlow::flow(){
	
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z0));

	su3_mat temp;
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
			temp = 0.25 * m_epsilon * (*m_Z0)(i, mu);
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
		}
	}
	
	(*m_U).transfer_FieldValues();
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z1));
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
			temp = m_epsilon *(8.0/9.0* (*m_Z1)(i, mu) -17.0/36.0*(*m_Z0)(i, mu));
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
		}
	}
	
	(*m_U).transfer_FieldValues();
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z2));
	for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
		for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
			temp = m_epsilon * (3.0 / 4.0 * (*m_Z2)(i, mu) - 8.0 / 9.0 * (*m_Z1)(i, mu)+ 17.0 / 36.0 * (*m_Z0)(i, mu));
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
		}
	}

	m_flowTime += m_epsilon;
	m_NrSteps++;
	(*m_U).transfer_FieldValues();

	if (m_NrSteps % m_FlowStepsPerMeasurement == 0) {
		MakeMeasurements();
	}

}
double GradientFlow::GetFlowTime() const{
	return m_flowTime;
}

void GradientFlow::Include_TopCharge(TopologicalCharge& TopCharge){
	(*m_U).transfer_FieldValues();
	m_topologicalCharge = &TopCharge;
}
void GradientFlow::Include_EnergyDensity(EnergyDensity& Edensity) {
	(*m_U).transfer_FieldValues();
	m_EnergyDensity = &Edensity;
}

void GradientFlow::MakeMeasurements(){
	if (m_topologicalCharge != nullptr) {
		m_topologicalCharge->calculate(m_flowTime);
	}
	if (m_EnergyDensity != nullptr) {
		m_EnergyDensity->calculate(m_flowTime);
	}
}
std::string GradientFlow::getupdateMethod() const {
	return m_updateMethod;
}