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
double Wilson::calculate_LocalActionChange(SU3_field& U_old, SU3_field& U_new, int i, int mu) {
	return -m_beta / 3.0 * ((U_new(i, mu) - U_old(i, mu)) * U_new.staple(i, mu)).ReTr();
}
void Wilson::calculate_LocalForce(SU3_field& field, SU3_field& F,int i, int mu) {
	su3_mat temp;

	su3_mat identity;
	identity.setToIdentity();
	temp = field(i, mu) * field.staple(i, mu);
	F(i, mu) = temp - temp.dagger();
	//F(i, mu) = field(i, mu) * Staple - (field(i, mu) * Staple).dagger();
	F(i, mu) = (m_beta / (6.0)) * ((0.5 * F(i, mu)) - (1.0 / 6.0) * (F(i, mu).Tr()) * identity).timesMinusI();
	if (!IsHermTrLess(F(i, mu))) {
		std::cout << "This F is not Hermitian and traceless!\n";
	}
}
void Wilson::calculate_Force(SU3_field& field, SU3_field& F){
	su3_mat identity;
	su3_mat temp;
	identity.setToIdentity();
	//RIPE FOR MULTI-THREADING
	for (int i = field.Responsible_Start(); i < field.Responsible_Stop(); i++) {
		for (int mu = 0; mu < field.getNrExtDOF(); mu++) {
			//Equation (8.42) in QCD on the lattice by Lang and Gattringer
			//Staple = field.staple(i, mu);

			//https://arxiv.org/pdf/1808.02281.pdf
			temp = field(i, mu) * field.staple(i, mu);
			F(i, mu) = temp - temp.dagger();
			//F(i, mu) = field(i, mu) * Staple - (field(i, mu) * Staple).dagger();
			F(i, mu) = (m_beta / (6.0)) * ((0.5 * F(i, mu)) - (1.0 / 6.0) * (F(i, mu).Tr()) * identity).timesMinusI();

#ifdef _DEBUG
			if (!IsHermTrLess(F(i, mu))) {
				std::cout << "This F is not Hermitian and traceless!\n";
			}
#endif // DEBUG
		}
	}
}

void Wilson::calculate_FlowGradient(SU3_field& field, SU3_field& F) {
	su3_mat identity;
	identity.setToIdentity();
	su3_mat temp;
	//RIPE FOR MULTI-THREADING
	for (int mu = 0; mu < field.getNrExtDOF(); mu++) {
		for (int i = field.Responsible_Start(); i < field.Responsible_Stop(); i++) {
		
			//Equation (8.42) in QCD on the lattice by Lang and Gattringer
			//Staple = field.staple(i, mu);

			//https://arxiv.org/pdf/1808.02281.pdf
			temp = field(i, mu) * field.staple(i, mu);

			F(i, mu) =  0.5*(temp-temp.dagger()) -(1.0 / 6.0) * (temp-temp.dagger()).Tr() * identity;
#ifdef _DEBUG
			if (!IsAntiHermTrLess(F(i, mu))) {
				std::cout << "This F is not Hermitian and traceless!\n";
			}
#endif // DEBUG
		}
	}
}




