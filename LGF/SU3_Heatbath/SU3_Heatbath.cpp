#include "SU3_Heatbath.h"
#include <assert.h>


SU3_Heatbath::SU3_Heatbath(SU3_field& U,double beta, int OR_per_HB) {
	m_updateMethod = "Heatbath_" + std::to_string(OR_per_HB) + "_ORperHB";
	m_OR_per_HB = OR_per_HB;
	m_beta = beta;
	//m_rand = new RNG_field(U.getLatticePtr(), 1);
	m_rand = new Random();
	m_U = &U;
	(*m_U).transfer_FieldValues();
}


void SU3_Heatbath::update(int n_iter) {
	for (int i = 0; i < n_iter; i++) {
		Cabibbo_Marinari((*m_U));
		OverRelaxation((*m_U));
		if (mpiWrapper::id() == 0) {
			std::cout << "step: " << i << std::endl;
		}
	}
}

//void SU3_Heatbath::update(int n_iter) {
//	for (int i = 0; i < n_iter; i++) {
//		su3_mat W;
//		C_double TSR[4];
//		C_double temp;
//
//		for (int parity = 0; parity <= 1; parity++) {
//			for (int mu = 0; mu < U.getNrExtDOF(); mu++) {
//				for (int i = U.Responsible_Start(parity); i < U.Responsible_Stop(parity); i++) {
//					for (int SubGroupNr = 0; SubGroupNr < 3; SubGroupNr++) {
//						int m, n;
//						n = SubGroupNr + 1;
//						m = n / 3;
//						n -= m;
//						//R*U
//						W = U(i, mu) * U.staple(i, mu);
//
//						SU2_Heatbath(TSR, W, i, SubGroupNr);
//						for (int col = 0; col < 3; col++) {
//							temp = TSR[0] * U(i, mu)(m, col) + TSR[1] * U(i, mu)(n, col);
//							U(i, mu)(n, col) = TSR[2] * U(i, mu)(m, col) + TSR[3] * U(i, mu)(n, col);
//							U(i, mu)(m, col) = temp;
//						}
//
//					}
//				}
//				U.transfer_FieldValues();
//			}
//		}
//	}
//}
void SU3_Heatbath::OverRelaxation(SU3_field& U) {

	su3_mat V;
	double a;
	double r[4];
	double temp_r[4];
	su3_mat UA_su3;
	C_double su2_sub[4];
	C_double temp;
	for (int step = 0; step < m_OR_per_HB; step++) {
		for (int parity = 0; parity <= 1; parity++) {
			for (int mu = 0; mu < U.getNrExtDOF(); mu++) {
				for (int i = U.Responsible_Start(parity); i < U.Responsible_Stop(parity); i++) {
					for (int SubGroupNr = 0; SubGroupNr < 3; SubGroupNr++) {
						int m, n;
						n = SubGroupNr + 1;
						m = n / 3;
						n -= m;
						UA_su3 = U(i, mu) * U.staple(i, mu);
						GetSU2submatrix_O4_rep(UA_su3, SubGroupNr, r);
						a = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
						if (a < 0.00001) {
							temp_r[0] = 1;
							temp_r[1] = 0;
							temp_r[2] = 0;
							temp_r[3] = 0;
						}
						else {
							double oneOvera = (1.0 / a);
							temp_r[0] = oneOvera * r[0];
							temp_r[1] = -oneOvera * r[1];
							temp_r[2] = -oneOvera * r[2];
							temp_r[3] = -oneOvera * r[3];
						}
						r[0] = temp_r[0] * temp_r[0] - temp_r[1] * temp_r[1] - temp_r[2] * temp_r[2] - temp_r[3] * temp_r[3];
						temp_r[0] *= 2;
						r[1] = temp_r[0] * temp_r[1];
						r[2] = temp_r[0] * temp_r[2];
						r[3] = temp_r[0] * temp_r[3];

						su2_sub[0] = C_double(r[0], r[3]); su2_sub[1] = C_double(r[2], r[1]);
						su2_sub[2] = C_double(-r[2], r[1]); su2_sub[3] = C_double(r[0], -r[3]);
						for (int col = 0; col < 3; col++) {
							temp = su2_sub[0] * U(i, mu)(m, col) + su2_sub[1] * U(i, mu)(n, col);
							U(i, mu)(n, col) = su2_sub[2] * U(i, mu)(m, col) + su2_sub[3] * U(i, mu)(n, col);
							U(i, mu)(m, col) = temp;
						}
					}
				}
				U.transfer_FieldValues();
			}
		}
	}
}
void SU3_Heatbath::Cabibbo_Marinari(SU3_field& U) {
	su3_mat W;
	C_double TSR[4];
	C_double temp;

	for (int parity = 0; parity <= 1; parity++) {
		for (int mu = 0; mu < U.getNrExtDOF(); mu++) {
			for (int i = U.Responsible_Start(parity); i < U.Responsible_Stop(parity); i++) {
				for (int SubGroupNr = 0; SubGroupNr < 3; SubGroupNr++) {
					int m, n;
					n = SubGroupNr + 1;
					m = n / 3;
					n -= m;
					//R*U
					W = U(i, mu) * U.staple(i, mu);

					

					SU2_Heatbath(TSR, W, i, SubGroupNr);
					for (int col = 0; col < 3; col++) {
						temp = TSR[0] * U(i, mu)(m, col) + TSR[1] * U(i, mu)(n, col);
						U(i, mu)(n, col) = TSR[2] * U(i, mu)(m, col) + TSR[3] * U(i, mu)(n, col);
						U(i, mu)(m, col) = temp;
					}
					////S*(R*U)
					//W = U(i, mu) * U.staple(i, mu);

					//TSR[0] = W(0, 0); TSR[1] = W(0, 2);
					//TSR[2] = W(2, 0); TSR[3] = W(2, 2);

					//heatbath_SU2(TSR, W, i, SubGroupNr);
					//for (int col = 0; col < 3; col++) {
					//	temp = TSR[0] * U(i, mu)(0, col) + TSR[1] * U(i, mu)(2, col);
					//	U(i, mu)(2, col) = TSR[2] * U(i, mu)(0, col) + TSR[3] * U(i, mu)(2, col);
					//	U(i, mu)(0, col) = temp;
					//}

					////T*(S*R*U)
					//W = U(i, mu) * U.staple(i, mu);

					//TSR[0] = W(1, 1); TSR[1] = W(1, 2);
					//TSR[2] = W(2, 1); TSR[3] = W(2, 2);

					//heatbath_SU2(TSR, W, i, SubGroupNr);
					//for (int col = 0; col < 3; col++) {
					//	temp = TSR[0] * U(i, mu)(1, col) + TSR[1] * U(i, mu)(2, col);
					//	U(i, mu)(2, col) = TSR[2] * U(i, mu)(1, col) + TSR[3] * U(i, mu)(2, col);
					//	U(i, mu)(1, col) = temp;
					//}

				}
			}
			U.transfer_FieldValues();
		}
	}
}

void SU3_Heatbath::SU2_Heatbath(C_double* TSR,su3_mat& W, int internIdx, int SubGroupNr) {
	static const double TwoPI = 8.0 * atan(1.0);
	double su2_O4[4];
	GetSU2submatrix_O4_rep(W, SubGroupNr, su2_O4);
	double X, X_prime, C, A;
	double abs_val_x;
	double x[4];
	C_double su2_sub[4];
	C_double V[4];
	double delta_bar, phi, sin_theta, cos_theta;
	double norm = sqrt(su2_O4[0] * su2_O4[0] + su2_O4[1] * su2_O4[1] + su2_O4[2] * su2_O4[2] + su2_O4[3] * su2_O4[3]);
    double alpha = (norm * m_beta/3.0);
	su2_sub[0] = C_double(su2_O4[0] / norm, -su2_O4[3] / norm);	su2_sub[1] = C_double(-su2_O4[2] / norm, -su2_O4[1] / norm);
	su2_sub[2] = C_double(su2_O4[2] / norm, -su2_O4[1] / norm);	su2_sub[3] = C_double(su2_O4[0] / norm, su2_O4[3] / norm);

	//See Improved heatbath method for Monte Carlo calculations in lattice gauge theories
	// https://www-sciencedirect-com.ezproxy.uio.no/science/article/pii/0370269385916326
	// p.6 (p.398)
	//do {
	//	X = -log(1.0-(*m_rand)(internIdx, 0).Uniform_Double());
	//	X_prime = -log(1.0-(*m_rand)(internIdx, 0).Uniform_Double());
	//	C = pow(cos(TwoPI * (*m_rand)(internIdx, 0).Uniform_Double()),2);
	//	A = X * C;
	//	delta_bar = (X_prime + A)/alpha;
	//} while (pow((*m_rand)(internIdx, 0).Uniform_Double(),2) > (1.0 - 0.5 * delta_bar));
 //    x[0] = 1.0 - delta_bar;
	// abs_val_x = sqrt(1.0 - x[0] * x[0]);
	// //See Gattringer and Lang p. 87
 //    cos_theta = 2.0 * (*m_rand)(internIdx, 0).Uniform_Double() - 1.0;
	// phi = TwoPI * (*m_rand)(internIdx, 0).Uniform_Double();
	do {
		X = -log(1.0 - (*m_rand).Uniform_Double());
		X_prime = -log(1.0 - (*m_rand).Uniform_Double());
		C = pow(cos(TwoPI * (*m_rand).Uniform_Double()), 2);
		A = X * C;
		delta_bar = (X_prime + A) / alpha;
	} while (pow((*m_rand).Uniform_Double(), 2) > (1.0 - 0.5 * delta_bar));
	x[0] = 1.0 - delta_bar;
	abs_val_x = sqrt(1.0 - x[0] * x[0]);
	//See Gattringer and Lang p. 87
	cos_theta = 2.0 * (*m_rand).Uniform_Double() - 1.0;
	phi = TwoPI * (*m_rand).Uniform_Double();
     
	sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	//Construct the SU(2) matrix V
	x[1] = abs_val_x * sin_theta * cos(phi);
	x[2] = abs_val_x * sin_theta * sin(phi);
	x[3] = abs_val_x * cos_theta;
	V[0] = C_double(x[0], x[3]);	V[1] = C_double(x[2], x[1]);
	V[2] = C_double(-x[2], x[1]);	V[3] = C_double(x[0], -x[3]);
	//update the subgroup by multiplication with V
	TSR[0] = V[0] * su2_sub[0] + V[1] * su2_sub[2];
	TSR[1] = V[0] * su2_sub[1] + V[1] * su2_sub[3];
	TSR[2] = V[2] * su2_sub[0] + V[3] * su2_sub[2];
	TSR[3] = V[2] * su2_sub[1] + V[3] * su2_sub[3];
}

std::string SU3_Heatbath::getupdateMethod() const{
	return m_updateMethod;
}

