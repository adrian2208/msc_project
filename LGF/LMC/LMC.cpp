#include "LMC.h"


LMC::LMC(SU3_field& U, double epsilon) {
	m_NrAccepted = 0;
	m_NrstepsTaken = 0;
	m_epsilon = epsilon;
	m_U = &U;
	exp_val_Plaquette = 0.0;
	exp_val_expdeltaH = 0.0;

}



void LMC::update() {
	double beta = 6.1;
	int N = m_U->getLatticePtr().m_responsible_Volume;
	su3_mat U_init;
	Random rng;
	su3_mat randomSU3;
	su3_mat A;
	for (int mu = 0; mu < m_U->getNrExtDOF(); mu++) {
	for (int i = 0; i < m_U->getLatticePtr().m_responsible_Volume; i++) {
		
			U_init = (*m_U)(i, mu);
			A = m_U->staple(i, mu);
			
			for (int hit = 0; hit < 20; hit++) {
				rng.rnd_su3_alg(randomSU3, m_epsilon);
				(*m_U)(i, mu) = randomSU3 * (*m_U)(i, mu);
				double expdeltaS = exp(-(beta / (double)N) * (((*m_U)(i, mu) - U_init) * A).ReTr());
				//std::cout << "expdeltaS = " << expdeltaS << "\n";
				if (rng.Uniform_Double() < expdeltaS) {
					m_NrAccepted++;
					exp_val_Plaquette += (*m_U).Avg_Plaquette();
					
				}
				else {
					(*m_U)(i, mu) = U_init;
				}
				m_NrstepsTaken++;
			}
			//std::cout << "Acceptance rate: " << (double)m_NrstepsTaken / ((double)m_NrAccepted) << "\n";
			//std::cout << "exp_val_Plaquette: " << exp_val_Plaquette / ((double)m_NrAccepted) << "\n\n";

		}
		
	}
	
}

double LMC::acceptanceRate() {
	return ((double)m_NrAccepted) / ((double)m_NrstepsTaken);
}