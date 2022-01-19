#pragma once
#include <random>
#include "SU3_mat.h"
#include <iostream>
class Random {

public:
	Random(int seed) {
		m_rng = std::mt19937_64(seed);
}
	Random() {
		std::random_device seed;
		m_rng = std::mt19937_64(seed());
	}
	double Uniform_Double() {
		std::uniform_real_distribution<double> out(0, 1);
		return out(m_rng);
	}
	double Gaussian_Double(double mean, double std) {
		std::normal_distribution<double> out(mean,std);
		return out(m_rng);
	}

	void rnd_su3_alg(su3_mat& mat, double epsilon) {
		
		//check whether its ok to use the same distribution many times...
		std::uniform_real_distribution<double> distribution(-0.5, 0.5);
		double norm, r0, r1, r2, r3, x0, x1, x2, x3;
		C_double u11[3], u12[3], u21[3], u22[3];
		for (int i = 0; i < 3; i++) {
			r0 = distribution(m_rng);
			r1 = distribution(m_rng);
			r2 = distribution(m_rng);
			r3 = distribution(m_rng);
			norm = sqrt(r1 * r1 + r2 * r2 + r3 * r3);

			x0 = sign(r0) * sqrt(1 - epsilon * epsilon);
			x1 = epsilon * r1 / norm;
			x2 = epsilon * r2 / norm;
			x3 = epsilon * r3 / norm;

			u11[i].R() = x0;
			u11[i].I() = x3;

			u12[i].R() = x2;
			u12[i].I() = x1;

			u21[i].R() = -x2;
			u21[i].I() = x1;

			u22[i].R() = x0;
			u22[i].I() = -x3;
		}
		mat[0] = u11[0] * u11[1];
		mat[1] = u11[0] * u12[1] * u21[2] + u12[0] * u11[2];
		mat[2] = u11[0] * u12[1] * u22[2] + u12[0] * u12[2];
		mat[3] = u21[0] * u11[1];
		mat[4] = u21[0] * u12[1] * u21[2] + u22[0] * u11[2];
		mat[5] = u21[0] * u12[1] * u22[2] + u22[0] * u12[2];
		mat[6] = u21[1];
		mat[7] = u22[1] * u21[2];
		mat[8] = u22[1] * u22[2];
		//std::cout << "r0:" << r0 << " " << "r1:" << r1 << " " << "r2:" << r2 << " " << "r3:" << r3 << "\n";
	}
private:
	std::mt19937_64 m_rng;

	inline int sign(double val) {
		return (val > 0) - (val < 0);
	}

};