#include "LHMC.h"


LHMC::LHMC(SU3_field& U, Wilson& GaugeAction, double epsilon) {
	m_NrAccepted = 0;
	m_NrstepsTaken = 0;
	m_epsilon = epsilon;
	m_NrleapfrogSteps = (int)(1.0 / abs(m_epsilon));
	m_GaugeAction = &GaugeAction;
	m_U = &U;
	m_P = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_F = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_U_init = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_F_init = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_rand = new RNG_field(U.getLatticePtr(), U.getNrExtDOF());
	m_K_init = new Double_field(U.getLatticePtr(),U.getNrExtDOF());
	m_localEpsilon = new Double_field(U.getLatticePtr(), U.getNrExtDOF());
	//m_S_init = new Double_field(U.getLatticePtr());
	(*m_U).transfer_FieldValues();
	(*m_U_init) = (*m_U);
	//m_GaugeAction->calculate_Force((*m_U), (*m_F));
	//for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
	//	for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
	//		(*m_S_init)(i, mu) = m_GaugeAction->calculate_LocalAction((*m_U),i,mu);
	//	}
	//}
	avg_Plaquette = 0.0;
#ifdef _DEBUG

	exp_val_expdeltaH = 0.0;
#endif // _DEBUG

}

double LHMC::calculate_KineticEnergy(SU3_field& P,int i, int mu) {
	return (P(i, mu) * P(i, mu)).ReTr();;
}

/// <summary>
/// Leapfrogs according to the procedure described in eq (8.38) in QCD on the lattice by Gattringer and Lang
/// </summary>
void LHMC::leapfrog(SU3_field& U, SU3_field& P, SU3_field& F, Wilson& action, int parity,int i,int mu) {

	su3_mat temp;
	//initial half-step
	m_GaugeAction->calculate_LocalForce(U, F, i, mu);
	P(i, mu) = P(i, mu) - (m_epsilon / 2.0) * F(i, mu);

	//n-2 full steps
	for (int k = 1; k < m_NrleapfrogSteps; k++) {
		//update U
		temp = P(i, mu).timesI() * m_epsilon;
		U(i, mu) = HermTrLessExp(temp) * U(i, mu);
		//update F
		//U.transfer_FieldValues();
		//action.calculate_Force(U, F);
		
		m_GaugeAction->calculate_LocalForce(U, F, i, mu);
		
		//update P
		P(i, mu) = P(i, mu) - m_epsilon * F(i, mu);
	}
	//final half-step
	//update U
	temp = P(i, mu).timesI() * m_epsilon;
	U(i, mu) = HermTrLessExp(temp) * U(i, mu);

	//update F
	//U.transfer_FieldValues();
	//action.calculate_Force(U, F);
	
	m_GaugeAction->calculate_LocalForce(U, F, i, mu);

	//update P
	P(i, mu) = P(i, mu) - (m_epsilon / 2.0) * F(i, mu);
}
void LHMC::sweep() {
	update(0);
	update(1);
	avg_Plaquette = (*m_U).Avg_Plaquette();
	std::cout << "exp_val_Plaquette: " << avg_Plaquette << "\n\n";
}
void LHMC::update(int parity) {
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	for (int mu = 0; mu < m_U->getNrExtDOF(); mu++) {
		(*m_U).transfer_FieldValues();
		for (int i = (*m_U).Responsible_Start(parity); i < (*m_U).Responsible_Stop(parity); i++) {
		
			(*m_U_init)(i, mu) = (*m_U)(i, mu);
			

			for (int j = 0; j < 8; j++) {
				LinCombGen = LinCombGen + rng.Gaussian_Double(0.0, 1.0) * generators(j);
			}
			(*m_P)(i, mu) = LinCombGen;
			LinCombGen.setToZeros();

			(*m_K_init)(i,mu) = calculate_KineticEnergy((*m_P),i,mu);

			//m_GaugeAction->calculate_LocalForce((*m_U), (*m_F), i, mu);
			//(*m_F_init)(i, mu) = (*m_F)(i, mu);
			leapfrog((*m_U), (*m_P), (*m_F), (*m_GaugeAction), parity,i,mu);

			double DeltaS = m_GaugeAction->calculate_LocalActionChange((*m_U_init), (*m_U), i, mu);//new-old
			double DeltaKE = calculate_KineticEnergy((*m_P), i, mu) - (*m_K_init)(i, mu);
			double DeltaH = DeltaS + DeltaKE;

			double Uni_double = (*m_rand)(i, mu).Uniform_Double();
			double exp_minDeltaH = exp(-DeltaH);
			//std::cout << "exp_minDeltaH = " << exp_minDeltaH << ", DeltaS = " << DeltaS << ", DeltaKE = " << DeltaKE << "\n";
			if (Uni_double < exp_minDeltaH) {
				m_NrAccepted++;
				//std::cout << "Accepted!\n";
			}
			else {
				(*m_U)(i, mu) = (*m_U_init)(i, mu);
				//std::cout << "Not Accepted!\n";
			}

			m_NrstepsTaken++;

		}
	}
	
}
	

double LHMC::acceptanceRate() {
	return ((double)m_NrAccepted) / ((double)m_NrstepsTaken);
}
