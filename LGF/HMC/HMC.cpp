#include "HMC.h"


HMC::HMC(SU3_field& U, Action& GaugeAction, double epsilon){
	m_NrAccepted = 0;
	m_NrstepsTaken = 0;
	m_epsilon = epsilon;
	m_NrleapfrogSteps = (int)1.0 / m_epsilon;
	m_GaugeAction = &GaugeAction;
	m_U = &U;
	m_P = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_F = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_U_init = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_F_init = new SU3_field(U.getLatticePtr(), U.getNrExtDOF());
	m_rand = new RNG_field(U.getLatticePtr());
	m_GaugeAction->calculate_Force((*m_U), (*m_F));
	S_init = m_GaugeAction->calculate_Action((*m_U));
}

double HMC::calculate_KineticEnergy(SU3_field& P) {
	double localSum = 0.0;
	double totalSum = 0.0;
	for (int i = 0; i < P.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < P.getNrExtDOF(); mu++) {
			localSum += (P(i, mu) * P(i, mu)).ReTr();
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
				temp = m_epsilon * P(i, mu);// .timesI();

#ifdef _DEBUG
				if (!IsHermTrLess(P(i, mu))) {
					std::cout << "^  P\n";
				}
				if (!IsHermTrLess(temp)) {
					std::cout << "^  temp\n";
				}
#endif // _DEBUG




				U(i, mu) = HermTrLessExp(temp) * U(i, mu);
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
			temp = m_epsilon * P(i, mu).timesI();
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
	m_U_init = m_U;
	m_F_init = m_F;
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	//double eps = 0.0001;//DELETE ME IF BELOW DOESNT WORK
	for (int i = 0; i < m_P->getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
			for (int j = 0; j < 8; j++) {
				LinCombGen = LinCombGen + rng.Uniform_Double() * generators(j);
			}
			//(*m_rand)(i, 0).rnd_su3_alg((*m_P)(i, mu), eps);//MAY NOT WORK.. SEE WORKSHEET
			(*m_P)(i, mu) = LinCombGen;
			LinCombGen.setToZeros();
		}
	}
	double H_init = S_init+ calculate_KineticEnergy((*m_P));
	leapfrog((*m_U), (*m_P), (*m_F), (*m_GaugeAction));


#ifdef _DEBUG
	for (int i = 0; i < m_P->getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < m_P->getNrExtDOF(); mu++) {
			if (!IsHermTrLess((*m_P)(i, mu))) {
				std::cout << "^  P after leapfrog\n";
			}
			if (!IsHermTrLess((*m_F)(i, mu))) {
				std::cout << "^  F after leapfrog\n";
			}
		}
	}
#endif // _DEBUG

	double S_final = m_GaugeAction->calculate_Action((*m_U));
	double H_final = S_final + calculate_KineticEnergy((*m_P));
	if (rng.Uniform_Double() < exp(H_final)) {
		S_init = S_final;
		m_NrAccepted++;
	}
	else {
		m_U = m_U_init;
		m_F = m_F_init;
	}
	m_NrstepsTaken++;
}
