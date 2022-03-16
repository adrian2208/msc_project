#include "LHMC.h"
#include <assert.h>
#include <algorithm>


LHMC::LHMC(SU3_field& U, Wilson& GaugeAction, double epsilon, int NrLeaps) {
	std::string eps_str = std::to_string(epsilon);
	std::replace(eps_str.begin(), eps_str.end(), '.', '_');
	m_updateMethod = "LHMC_" + std::to_string(NrLeaps) + "_of_"+ eps_str;
	m_NrAccepted = 0;
	m_NrstepsTaken = 0;
	m_epsilon = epsilon;
	if (NrLeaps == 0) {
		m_NrleapfrogSteps = (int)(1.0 / abs(m_epsilon));
		assert(m_NrleapfrogSteps != 0);
	}
	else {
		m_NrleapfrogSteps = NrLeaps;
		assert(m_NrleapfrogSteps*epsilon < 2 && m_NrleapfrogSteps * epsilon >0.5);
	}
	m_GaugeAction = &GaugeAction;
	m_U = &U;
	m_rand = new Random();//new RNG_field(U.getLatticePtr(), 1);
	(*m_U).transfer_FieldValues();

	avg_Plaquette = 0.0;
#ifdef _DEBUG

	exp_val_expdeltaH = 0.0;
#endif // _DEBUG

}
LHMC::~LHMC() {
	delete m_rand;
}

double LHMC::calculate_KineticEnergy(su3_mat& P) {
	return (P * P).ReTr();;
}

/// <summary>
/// Leapfrogs according to the procedure described in eq (8.38) in QCD on the lattice by Gattringer and Lang
/// </summary>
void LHMC::leapfrog(SU3_field& U, su3_mat& P, su3_mat& F, Wilson& action, int parity,int i,int mu) {
	double beta = action.getBeta();
	su3_mat temp;
	su3_mat staple;
	su3_mat Ftemp;

	su3_mat identity;
	identity.setToIdentity();
	staple = U.staple(i, mu);
	//initial half-step
	//m_GaugeAction->calculate_LocalForce(U, F, i, mu);
	Ftemp = U(i, mu) * staple;
	F = Ftemp - Ftemp.dagger();
	F = (beta / (6.0)) * ((0.5 * F) - (1.0 / 6.0) * (F.Tr()) * identity).timesMinusI();
	
	P = P - (m_epsilon / 2.0) * F;

	//n-2 full steps
	for (int k = 1; k < m_NrleapfrogSteps; k++) {
		//update U
		temp = P.timesI() * m_epsilon;
		U(i, mu) = HermTrLessExp(temp) * U(i, mu);
		//update F
		
		//m_GaugeAction->calculate_LocalForce(U, F, i, mu);
		Ftemp = U(i, mu) * staple;
		F = Ftemp - Ftemp.dagger();
		F = (beta / (6.0)) * ((0.5 * F) - (1.0 / 6.0) * (F.Tr()) * identity).timesMinusI();
		
		//update P
		P = P - m_epsilon * F;
	}
	//final half-step
	//update U
	temp = P.timesI() * m_epsilon;
	U(i, mu) = HermTrLessExp(temp) * U(i, mu);

	//update F
	//m_GaugeAction->calculate_LocalForce(U, F, i, mu);
	Ftemp = U(i, mu) * staple;
	F = Ftemp - Ftemp.dagger();
	F = (beta / (6.0)) * ((0.5 * F) - (1.0 / 6.0) * (F.Tr()) * identity).timesMinusI();

	//update P
	P = P - (m_epsilon / 2.0) * F;
}
void LHMC::sweep() {
	update(0);
	MPI_Barrier(mpiWrapper::comm());
	update(1);
	MPI_Barrier(mpiWrapper::comm());
}
void LHMC::update(int parity) {

	//
	// THIS IS WRONG! MU -> EVEN -> ODD .... NOT EVEN -> MU ->NU ->... -> ODD -> ....
	// 
	//su3_mat P;
	//su3_mat F;
	//SU3_gen generators;
	//double KE_init;
	//su3_mat U_init;
	//for (int mu = 0; mu < m_U->getNrExtDOF(); mu++) {
	//	(*m_U).transfer_FieldValues();
	//	for (int i = (*m_U).Responsible_Start(parity); i < (*m_U).Responsible_Stop(parity); i++) {
	//	
	//		U_init = (*m_U)(i, mu);
	//		
	//		P.setToZeros();
	//		for (int j = 0; j < 8; j++) {
	//			P = P + (*m_rand).Gaussian_Double(0.0, 1.0) * generators(j);
	//		}

	//		KE_init = calculate_KineticEnergy(P);

	//		leapfrog((*m_U), P, F, (*m_GaugeAction), parity,i,mu);

	//		double DeltaS = m_GaugeAction->calculate_LocalActionChange(U_init, (*m_U), i, mu);//new-old
	//		double DeltaKE = calculate_KineticEnergy(P) - KE_init;
	//		double DeltaH = DeltaS + DeltaKE;

	//		double Uni_double = (*m_rand).Uniform_Double();
	//		double exp_minDeltaH = exp(-DeltaH);
	//		//std::cout << "exp_minDeltaH = " << exp_minDeltaH << ", DeltaS = " << DeltaS << ", DeltaKE = " << DeltaKE << "\n";
	//		if (Uni_double < exp_minDeltaH) {
	//			m_NrAccepted++;
	//			//std::cout << "Accepted!\n";
	//		}
	//		else {
	//			(*m_U)(i, mu) = U_init;
	//			//std::cout << "Not Accepted!\n";
	//		}

	//		m_NrstepsTaken++;

	//	}
	//}
	
}
	

double LHMC::acceptanceRate() {
	return ((double)m_NrAccepted) / ((double)m_NrstepsTaken);
}
std::string LHMC::getupdateMethod() const {
	return m_updateMethod;
}