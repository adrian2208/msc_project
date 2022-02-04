#include "SU3_Heatbath.h"
#include <assert.h>


SU3_Heatbath::SU3_Heatbath(SU3_field& U) {

	m_rand = new RNG_field(U.getLatticePtr(), 1);

}


void SU3_Heatbath::update(SU3_field& U, int n_iter, double beta) {

	double zeta;
	zeta = 1;
	int i, j, k, iter, mu, parity;

	su3_mat M;
	C_double a[4], tmpUik;

	for (iter = 0; iter < n_iter; iter++)
		for (parity = 0; parity < 2; parity++)
			for (mu = 0; mu < U.getNrExtDOF(); mu++) {
				for (int site = U.Responsible_Start(parity); site < U.Responsible_Stop(parity); site++) {
					for (i = 0; i < 3 - 1; i++)
						for (j = i + 1; j < 3; j++) {
							if (zeta == 1)		M = U(site, mu) * U.staple(site,mu);
							else if (mu == 0)	M = zeta * U(site, 0) * U.staple(site, 0);
							else				M = ((double)1.0 / zeta) * U(site, mu) * U.staple(site, mu);
							a[0] = M(i, i);
							a[1] = M(i, j);
							a[2] = M(j, i);
							a[3] = M(j, j);
							heatbath_SU2((*m_rand)(site,0), beta / 3.0, a);
							for (k = 0; k < 3; k++) {
								tmpUik = a[0] * U(site, mu)(i, k) + a[1] * U(site, mu)(j, k);
								U(site, mu)(j, k) = a[2] * U(site, mu)(i, k) + a[3] * U(site, mu)(j, k);
								U(site, mu)(i, k) = tmpUik;
							}
						}
				}
				U.transfer_FieldValues();
			}
}

void SU3_Heatbath::heatbath_SU2(Random& random, double beta_eff, C_double* a) {
	static const double PI = 4.0 * atan(1.0);
	double e0, e1, e2, e3, dk, p0;
	double r1, r2, r3, r4;
	double a0, a1, a2, a3;
	C_double u0, u1, u2, u3;
	C_double v0, v1, v2, v3;
	double delta, phi, sin_alpha, sin_theta, cos_theta;
    e0 = a[0].R() + a[3].R();
    e1 = a[1].I() + a[2].I();
    e2 = a[1].R() - a[2].R();
    e3 = a[0].I() - a[3].I();
    dk = sqrt(e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3);
    p0 = (dk * beta_eff);
    u0 = C_double(e0 / dk, -e3 / dk);
    u2 = C_double(e2 / dk, -e1 / dk);
    u1 = C_double(-e2 / dk, -e1 / dk);
    u3 = C_double(e0 / dk, e3 / dk);

	do {
       do; while ((r1 = random.Uniform_Double()) < 0.0001);
       r1 = -log(r1) / p0;
       do; while ((r2 = random.Uniform_Double()) < 0.0001);
       r2 = -log(r2) / p0;
       r3 = cos(2.0 * PI * random.Uniform_Double());
       r3 = r3 * r3;
       delta = r2 + r1 * r3;
       r4 = random.Uniform_Double();
	} while (r4 * r4 > (1.0 - 0.5 * delta));
     a0 = 1.0 - delta;
     cos_theta = 2.0 * random.Uniform_Double() - 1.0;
     sin_theta = sqrt(1.0 - cos_theta * cos_theta);
     sin_alpha = sqrt(1 - a0 * a0);
     phi = 2.0 * PI * random.Uniform_Double();
     a1 = sin_alpha * sin_theta * cos(phi);
     a2 = sin_alpha * sin_theta * sin(phi);
     a3 = sin_alpha * cos_theta;
     v0 = C_double(a0, a3);
     v1 = C_double(a2, a1);
     v2 = C_double(-a2, a1);
     v3 = C_double(a0, -a3);
     a[0] = v0 * u0 + v1 * u2;
     a[1] = v0 * u1 + v1 * u3;
     a[2] = v2 * u0 + v3 * u2;
     a[3] = v2 * u1 + v3 * u3;
}