#include "GradientFlow.h"

GradientFlow::GradientFlow(Action& GaugeAction, SU3_field& U, double epsilon){
	m_flowTime = 0.0;
	m_GaugeAction = &GaugeAction;
	m_epsilon = epsilon;
	m_U = &U;
	m_Z0 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_Z1 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_Z2 = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
}

void GradientFlow::flow(){
	(*m_U).transfer_FieldValues();
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z0));



	su3_mat temp;
	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
		for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
			temp = -0.25 * m_epsilon * (*m_Z0)(i, mu);
#ifdef _DEBUG
			if (!isSpecialUnitary((*m_U)(i, mu))) {
				std::cout << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS ^ \n";
			}
#endif//debug
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
			#ifdef _DEBUG
			if (!isSpecialUnitary((*m_U)(i, mu))) {
				std::cout << "W_1 ^ \n";
			}
			#endif//debug
		}
	}
	
	(*m_U).transfer_FieldValues();
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z1));

	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
		for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
			temp = -m_epsilon *(8.0/9.0* (*m_Z1)(i, mu) -17.0/36.0*(*m_Z0)(i, mu));
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
			#ifdef _DEBUG
			if (!isSpecialUnitary((*m_U)(i, mu))) {
				std::cout << "W_2 ^ \n";
			}
			#endif//debug
		}
	}
	
	(*m_U).transfer_FieldValues();
	(*m_GaugeAction).calculate_FlowGradient((*m_U), (*m_Z2));

	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
		for (int mu = 0; mu < (*m_U).getNrExtDOF(); mu++) {
			temp = -m_epsilon * (3.0 / 4.0 * (*m_Z2)(i, mu) - 8.0 / 9.0 * (*m_Z1)(i, mu)+ 17.0 / 36.0 * (*m_Z0)(i, mu));
			(*m_U)(i, mu) = HermTrLessExp(temp) * (*m_U)(i, mu);
			#ifdef _DEBUG
			if (!isSpecialUnitary((*m_U)(i, mu))) {
				std::cout << "W_3 ^ \n";
			}
			#endif//debug
		}
	}

	m_flowTime += m_epsilon;
	MakeMeasurements();
}

void GradientFlow::Include_TopCharge(TopologicalCharge& TopCharge){
	(*m_U).transfer_FieldValues();
	m_topologicalCharge = &TopCharge;
}

void GradientFlow::MakeMeasurements(){
	if (m_topologicalCharge != nullptr) {
		m_topologicalCharge->calculate();
	}
}
