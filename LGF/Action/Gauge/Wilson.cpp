#include "Wilson.h"

Wilson::Wilson(double beta) : Action() {
	m_beta = beta;
	//m_Coeff = field.getLatticePtr().m_totalVolume * field.getNrExtDOF() * (field.getNrExtDOF() - 1) * beta / 2.0;
}

double Wilson::calculate_Action(SU3_field& field){
	return field.getLatticePtr().m_totalVolume * field.getNrExtDOF() * (field.getNrExtDOF() - 1) * m_beta / 2.0 *1.0-m_beta*field.total_PlaquetteSum()/2.0;
}

void Wilson::calculate_Force(SU3_field& field, SU3_field& F){
	su3_mat Staple;
	for (int i = 0; i < field.getLatticePtr().m_responsible_Volume; i++) {
		for (int mu = 0; mu < field.getNrExtDOF(); mu++) {
			//Equation (8.42) in QCD on the lattice by Lang and Gattringer
			Staple = field.staple(i, mu);
			F(i, mu) = m_beta/12.0*(field(i, mu) *Staple -(field(i, mu)* Staple).dagger()).timesMinusI();
		}
	}
}
