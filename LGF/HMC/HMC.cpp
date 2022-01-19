#include "HMC.h"


HMC::HMC(SU3_field& U, Action& GaugeAction, double epsilon){
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
	m_rand = new RNG_field(U.getLatticePtr());
	m_GaugeAction->calculate_Force((*m_U), (*m_F));
	S_init = m_GaugeAction->calculate_Action((*m_U));


	exp_val_Plaquette = 0.0;
#ifdef _DEBUG
	
	exp_val_expdeltaH = 0.0;
#endif // _DEBUG

}

double HMC::calculate_KineticEnergy(SU3_field& P) {
	double localSum = 0.0;
	double totalSum = 0.0;
	for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
			localSum +=(P(i, mu) * P(i, mu)).ReTr();
		}
	}
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	return totalSum;
}

/// <summary>
/// Leapfrogs according to the procedure described in eq (8.38) in QCD on the lattice by Gattringer and Lang
/// </summary>
void HMC::leapfrog(SU3_field& U, SU3_field& P, SU3_field& F,Action& action){

	su3_mat temp;
	//initial half-step
	for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
			P(i,mu) = P(i, mu) - (m_epsilon/2.0) * F(i, mu);
		}
	}


#ifdef _DEBUG
	for (int i = 0; i < m_P->getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
			if (!IsHermTrLess(P(i, mu))) {
				std::cout << "^  P while leapfrog\n";
			}
			if (!IsHermTrLess(F(i, mu))) {
				std::cout << "^  F while leapfrog\n";
			}
		}
	}

#endif // _DEBUG
	//n-2 full steps
	for (int k = 1; k < m_NrleapfrogSteps; k++) {
		//update U
		for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
			for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
				temp = P(i, mu).timesI()*m_epsilon;
				U(i, mu) = HermTrLessExp(temp) * U(i, mu);
#ifdef _DEBUG
				bool somethingWrong = false;
				if (!IsHermTrLess(P(i, mu))) {
					std::cout << "^  P\n";
					somethingWrong = true;
				}
				if (!IsAntiHermTrLess(temp)) {
					std::cout << "^  temp\n";
					somethingWrong = true;
				}
				if (!isSpecialUnitary(U(i, mu))) {
					std::cout << "^  U\n";
					somethingWrong = true;
				}
				if (somethingWrong) {
					std::cout << "Error on step nr. " << m_NrstepsTaken << "\n";
				}
#endif // _DEBUG
			}
		}
		//update F
		action.calculate_Force(U, F);
		//update P
		for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
			for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
				P(i, mu) = P(i, mu) - m_epsilon * F(i, mu);
			}
		}
	}
	//final half-step
	//update U
	for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
			temp = P(i, mu).timesI()*m_epsilon;
			U(i, mu) = HermTrLessExp(temp) * U(i, mu);
		}
	}

	//update F
	action.calculate_Force(U, F);
	//update P
	for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
			P(i, mu) = P(i, mu) - (m_epsilon/2.0) * F(i, mu);
		}
	}

}

void HMC::update(){
	(*m_U_init)=(*m_U);
	(*m_F_init) = (*m_F);
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	//double eps = 0.0001;//DELETE ME IF BELOW DOESNT WORK
	for (int i = 0; i < m_P->getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
			for (int j = 0; j < 8; j++) {
				LinCombGen = LinCombGen + rng.Gaussian_Double(0.0,1.0) * generators(j);
			}
			//(*m_rand)(i, 0).rnd_su3_alg((*m_P)(i, mu), eps);//MAY NOT WORK.. SEE WORKSHEET
			(*m_P)(i, mu) = LinCombGen;
			LinCombGen.setToZeros();
		}
	}
	double KE_init = calculate_KineticEnergy((*m_P));
	double H_init = S_init + KE_init;
	leapfrog((*m_U), (*m_P), (*m_F), (*m_GaugeAction));

#ifdef _DEBUG
	for (int i = 0; i < m_P->getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
			bool somethingWrong = false;
			if (!IsHermTrLess((*m_P)(i, mu))) {
				std::cout << "^  P\n";
				somethingWrong = true;
			}
			if (!IsHermTrLess((*m_F)(i, mu))) {
				std::cout << "^  F\n";
				somethingWrong = true;
			}
			if (!isSpecialUnitary((*m_U)(i, mu))) {
				//std::cout << "^  U\n";
				somethingWrong = true;
			}
			if (somethingWrong) {
				//std::cout << "Error after leapfrog on step nr. " << m_NrstepsTaken << "\n";
			}
		}
	}
#endif // _DEBUG

	double S_final = m_GaugeAction->calculate_Action((*m_U));
	double KE = calculate_KineticEnergy((*m_P));
	double H_final = S_final + KE;

	std::cout << "S_final = " << S_final << ", S_initial = " << S_init << "\n";
	std::cout << "KE_final = " << KE << ", KE_initial = " << KE_init << "\n";
	std::cout << "H_final = " << H_final << ", H_initial = " << H_init << "\n\n\n";
#ifdef _DEBUG
	exp_val_expdeltaH+= exp(H_init - H_final);
	std::cout << "exp_val_expdeltaH: " << exp_val_expdeltaH / ((double)(m_NrstepsTaken+1)) << "\n";
#endif // _DEBUG
	
	if (rng.Uniform_Double() < exp(H_init- H_final)) {
		S_init = S_final;
		m_NrAccepted++;

		exp_val_Plaquette+= (*m_U).Avg_Plaquette();
		std::cout << "exp_val_Plaquette: " << exp_val_Plaquette / ((double)m_NrAccepted) << "\n\n";

	}
	else {
		(*m_U) = (*m_U_init);
		(*m_F) = (*m_F_init);
	}
	m_NrstepsTaken++;
}

double HMC::acceptanceRate(){
	return ((double)m_NrAccepted)/((double)m_NrstepsTaken);
}
