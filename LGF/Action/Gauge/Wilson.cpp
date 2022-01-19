#include "Wilson.h"
//#define DEBUG
Wilson::Wilson(double beta) : Action() {
	m_beta = beta;
	//m_Coeff = field.getLatticePtr().m_totalVolume * field.getNrExtDOF() * (field.getNrExtDOF() - 1) * beta / 2.0;
}

double Wilson::calculate_Action(SU3_field& field){
	//return  m_beta* field.getLatticePtr().m_totalVolume * field.getNrExtDOF() * (field.getNrExtDOF() - 1) -m_beta * field.total_PlaquetteSum() / 3.0;//+field.getLatticePtr().m_totalVolume * field.getNrExtDOF() * (field.getNrExtDOF() - 1) * m_beta / 2.0 * 1.0;
	return m_beta * field.total_PlaquetteSum();
}

void Wilson::calculate_Force(SU3_field& field, SU3_field& F){
	su3_mat Staple;
	su3_mat identity;
	identity.setToIdentity();
	//RIPE FOR MULTI-THREADING
	for (int i = 0; i < field.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < field.getNrExtDOF(); mu++) {
			//Equation (8.42) in QCD on the lattice by Lang and Gattringer
			Staple = field.staple(i, mu);
			//F(i, mu) = field(i, mu) *Staple.dagger();
			//F(i, mu) = m_beta / 6.0 * (F(i, mu).at()).timesI();
			
			//https://arxiv.org/pdf/1808.02281.pdf
			F(i, mu) = field(i, mu) * Staple - (field(i, mu) * Staple).dagger();
			F(i, mu) = (m_beta /(6.0)) * ((0.5 * F(i, mu)) - (1.0 / 6.0) * (F(i, mu).Tr()) * identity).timesMinusI();
#ifdef _DEBUG
			if (!IsHermTrLess(F(i, mu))) {
				std::cout << "This F is not Hermitian and traceless!\n";
			}
#endif // DEBUG
		}
	}
}



