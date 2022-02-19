#include "SU3_Heatbath.h"
#include <assert.h>


SU3_Heatbath::SU3_Heatbath(SU3_field& U,double beta, int OR_per_HB) {
	m_updateMethod = "Heatbath_" + std::to_string(OR_per_HB) + "_ORperHB";
	m_OR_per_HB = OR_per_HB;
	m_beta = beta;
	m_rand = new RNG_field(U.getLatticePtr(), 1);
	m_U = &U;
	(*m_U).transfer_FieldValues();
}


void SU3_Heatbath::update(int n_iter) {
	for (int i = 0; i < n_iter; i++) {
		Cabibbo_Marinari((*m_U));
		OverRelaxation((*m_U));
	}
}
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
				//R*U
				W = U(i, mu) * U.staple(i, mu);

				TSR[0] = W(0,0); TSR[1] = W(0,1);
				TSR[2] = W(1,0); TSR[3] = W(1,1);

				heatbath_SU2(TSR, i);
				for (int col = 0; col < 3; col++) {
					temp = TSR[0] * U(i, mu)(0, col) + TSR[1] * U(i, mu)(1, col);
					U(i, mu)(1, col) = TSR[2] * U(i, mu)(0, col) + TSR[3] * U(i, mu)(1, col);
					U(i, mu)(0, col) = temp;
				}
				//S*(R*U)
				W = U(i, mu) * U.staple(i, mu);

				TSR[0] = W(0, 0); TSR[1] = W(0, 2);
				TSR[2] = W(2, 0); TSR[3] = W(2, 2);

				heatbath_SU2(TSR, i);
				for (int col = 0; col < 3; col++) {
					temp = TSR[0] * U(i, mu)(0, col) + TSR[1] * U(i, mu)(2, col);
					U(i, mu)(2, col) = TSR[2] * U(i, mu)(0, col) + TSR[3] * U(i, mu)(2, col);
					U(i, mu)(0, col) = temp;
				}

				//T*(S*R*U)
				W = U(i, mu) * U.staple(i, mu);

				TSR[0] = W(1,1); TSR[1] = W(1,2);
				TSR[2] = W(2,1); TSR[3] = W(2,2);

				heatbath_SU2(TSR, i);
				for (int col = 0; col < 3; col++) {
					temp = TSR[0] * U(i, mu)(1, col) + TSR[1] * U(i, mu)(2, col);
					U(i, mu)(2, col) = TSR[2] * U(i, mu)(1, col) + TSR[3] * U(i, mu)(2, col);
					U(i, mu)(1, col) = temp;
				}


			}
			U.transfer_FieldValues();
		}
	}
}

void SU3_Heatbath::heatbath_SU2(C_double* TSR, int internIdx) {
	static const double PI = 4.0 * atan(1.0);
	double e0, e1, e2, e3, dk, p0;
	double r1, r2, r3, r4;
	double a0, a1, a2, a3;
	C_double u0, u1, u2, u3;
	C_double v0, v1, v2, v3;
	double delta, phi, sin_alpha, sin_theta, cos_theta;
    e0 = TSR[0].R() + TSR[3].R();
    e1 = TSR[1].I() + TSR[2].I();
    e2 = TSR[1].R() - TSR[2].R();
    e3 = TSR[0].I() - TSR[3].I();
    dk = sqrt(e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3);
    p0 = (dk * m_beta/3.0);
    u0 = C_double(e0 / dk, -e3 / dk);
    u2 = C_double(e2 / dk, -e1 / dk);
    u1 = C_double(-e2 / dk, -e1 / dk);
    u3 = C_double(e0 / dk, e3 / dk);

	do {
       do; while ((r1 = (*m_rand)(internIdx,0).Uniform_Double()) < 0.0001);
       r1 = -log(r1) / p0;
       do; while ((r2 = (*m_rand)(internIdx, 0).Uniform_Double()) < 0.0001);
       r2 = -log(r2) / p0;
       r3 = cos(2.0 * PI * (*m_rand)(internIdx, 0).Uniform_Double());
       r3 = r3 * r3;
       delta = r2 + r1 * r3;
       r4 = (*m_rand)(internIdx, 0).Uniform_Double();
	} while (r4 * r4 > (1.0 - 0.5 * delta));
     a0 = 1.0 - delta;
     cos_theta = 2.0 * (*m_rand)(internIdx, 0).Uniform_Double() - 1.0;
     sin_theta = sqrt(1.0 - cos_theta * cos_theta);
     sin_alpha = sqrt(1 - a0 * a0);
     phi = 2.0 * PI * (*m_rand)(internIdx, 0).Uniform_Double();
     a1 = sin_alpha * sin_theta * cos(phi);
     a2 = sin_alpha * sin_theta * sin(phi);
     a3 = sin_alpha * cos_theta;
     v0 = C_double(a0, a3);
     v1 = C_double(a2, a1);
     v2 = C_double(-a2, a1);
     v3 = C_double(a0, -a3);
	 TSR[0] = v0 * u0 + v1 * u2;
	 TSR[1] = v0 * u1 + v1 * u3;
	 TSR[2] = v2 * u0 + v3 * u2;
	 TSR[3] = v2 * u1 + v3 * u3;
}

std::string SU3_Heatbath::getupdateMethod() const{
	return m_updateMethod;
}
