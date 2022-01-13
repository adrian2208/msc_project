#include "SU3_mat.h"
#include <cmath>
using namespace std;
su3_mat::su3_mat(){
	allocate_matrix();
}
su3_mat::su3_mat(const su3_mat& matrix) {
	allocate_matrix();
	for (int i = 0; i < 9; i++) {
		mat[i] = matrix[i];
	}
}

su3_mat::~su3_mat(){
	delete[] mat;
}

inline void su3_mat::allocate_matrix(){
	mat = new C_double[9];
}

inline C_double& su3_mat::operator() (int row, int col) {
	return mat[3 * row + col];
}

inline C_double& su3_mat::operator[] (int i) const{
	return mat[i];
}
inline const su3_mat& su3_mat::operator=(const su3_mat& a) {
	for (int i = 0; i < 9; i++) {
		mat[i] = a[i];
	}
	return *this;
}

su3_mat su3_mat::dagger() const{
	su3_mat out;
	out[0] = mat[0].dagger();
	out[1] = mat[3].dagger();
	out[2] = mat[6].dagger();
	out[3] = mat[1].dagger();
	out[4] = mat[4].dagger();
	out[5] = mat[7].dagger();
	out[6] = mat[2].dagger();
	out[7] = mat[5].dagger();
	out[8] = mat[8].dagger();
	return out;
}
su3_mat su3_mat::dagger(){
	su3_mat out;
	out[0] = mat[0].dagger();
	out[1] = mat[3].dagger();
	out[2] = mat[6].dagger();
	out[3] = mat[1].dagger();
	out[4] = mat[4].dagger();
	out[5] = mat[7].dagger();
	out[6] = mat[2].dagger();
	out[7] = mat[5].dagger();
	out[8] = mat[8].dagger();
	return out;
}

su3_mat su3_mat::timesI() const{
	su3_mat out;
	out[0] = mat[0].timesI();
	out[1] = mat[3].timesI();
	out[2] = mat[6].timesI();
	out[3] = mat[1].timesI();
	out[4] = mat[4].timesI();
	out[5] = mat[7].timesI();
	out[6] = mat[2].timesI();
	out[7] = mat[5].timesI();
	out[8] = mat[8].timesI();
	return out;
}

su3_mat su3_mat::timesMinusI() const{
	su3_mat out;
	out[0] = mat[0].timesMinusI();
	out[1] = mat[3].timesMinusI();
	out[2] = mat[6].timesMinusI();
	out[3] = mat[1].timesMinusI();
	out[4] = mat[4].timesMinusI();
	out[5] = mat[7].timesMinusI();
	out[6] = mat[2].timesMinusI();
	out[7] = mat[5].timesMinusI();
	out[8] = mat[8].timesMinusI();
	return out;
}

C_double su3_mat::det() const{	
	return mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]) - mat[1] * (mat[3] * mat[8] - mat[5] * mat[6]) + mat[2] * (mat[3] * mat[7] - mat[4] * mat[6]);
}

C_double su3_mat::Tr() const{
	return mat[0]+mat[4]+mat[8];
}

double su3_mat::ReTr() const{
	return mat[0].R() + mat[4].R() + mat[8].R();
}

void su3_mat::setToIdentity(){
	mat[0].Re = mat[4].Re = mat[8].Re = 1.0;
	mat[0].Im = mat[4].Im = mat[8].Im = 0.0;
	mat[1] = mat[2] = mat[3] = mat[5] = mat[6] = mat[7] = C_double(0.0, 0.0);
}
void su3_mat::setToZeros() {
	mat[0] = mat[1] = mat[2] = mat[3] = mat[4] = mat[5] = mat[6] = mat[7] = mat[8] = C_double(0.0, 0.0);
}

C_double* su3_mat::getMemPointer(){
	return mat;
}

su3_mat operator*(const su3_mat& a, const su3_mat& b){
	su3_mat out;
	out[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
	out[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
	out[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];
	out[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
	out[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
	out[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];
	out[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
	out[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
	out[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
	return out;
}

su3_mat operator+(const su3_mat& a, const su3_mat& b){
	su3_mat out;
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	out[2] = a[2] + b[2];
	out[3] = a[3] + b[3];
	out[4] = a[4] + b[4];
	out[5] = a[5] + b[5];
	out[6] = a[6] + b[6];
	out[7] = a[7] + b[7];
	out[8] = a[8] + b[8];
	return out;
}
su3_mat operator-(const su3_mat& a, const su3_mat& b) {
	su3_mat out;
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	out[2] = a[2] - b[2];
	out[3] = a[3] - b[3];
	out[4] = a[4] - b[4];
	out[5] = a[5] - b[5];
	out[6] = a[6] - b[6];
	out[7] = a[7] - b[7];
	out[8] = a[8] - b[8];
	return out;
}

su3_mat operator*(const su3_mat& a, const int b){
	su3_mat out;
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	out[3] = a[3] * b;
	out[4] = a[4] * b;
	out[5] = a[5] * b;
	out[6] = a[6] * b;
	out[7] = a[7] * b;
	out[8] = a[8] * b;
	return out;
}

su3_mat operator*(const int b, const su3_mat& a){
	return a*b;
}

su3_mat operator*(const su3_mat& a, const double b){
	su3_mat out;
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	out[3] = a[3] * b;
	out[4] = a[4] * b;
	out[5] = a[5] * b;
	out[6] = a[6] * b;
	out[7] = a[7] * b;
	out[8] = a[8] * b;
	return out;
}

su3_mat operator*(const double b, const su3_mat& a){
	return a*b;
}

su3_mat operator*(const su3_mat& a, const C_double b){
	su3_mat out;
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	out[3] = a[3] * b;
	out[4] = a[4] * b;
	out[5] = a[5] * b;
	out[6] = a[6] * b;
	out[7] = a[7] * b;
	out[8] = a[8] * b;
	return out;
}

su3_mat operator*(const C_double b, const su3_mat& a){
	return a*b;
}

su3_mat HermTrLessExp(su3_mat& iP){

#ifdef _DEBUG
	bool Is_iP_appropriate = IsHermTrLess(iP, true);
	if (!Is_iP_appropriate) {
		std::cout << "Attempted exponentiation of unfit matrix" << "\n";
	}
#endif // _DEBUG

	su3_mat Q = iP.timesMinusI();
	su3_mat Q_squared = Q * Q;
	double c0 = (Q * Q_squared).ReTr() / 3.0;//eq(14)
	double c0_abs = abs(c0);//explanation below eq(34)
	double c1 = Q_squared.ReTr() / 2.0;//eq(15)
	double c1_sqrt = sqrt(c1);
	double c0_max = 2.0 * pow((c1 / 3.0), 1.5);//eq(17)
	double absOverMax = c0_abs / c0_max;
	if (absOverMax > 1.0 || absOverMax< 0.0) {
		std::cout << "ERROR: attempted acos(x) with invalid entry: x = "<<absOverMax << "\n";
	}
	double theta_div3 = acos(absOverMax) / 3.0;//eq(25)
	double u = sqrt(1.0 / 3.0)*c1_sqrt * cos(theta_div3);//eq(23)
	double u_squared = u * u;
	double u_cos = cos(u);
	double u_sin = sin(u);
	double w = c1_sqrt * sin(theta_div3);//eq(23)
	double w_squared = w * w;
	double w_cos = cos(w);
	double xi0;//between eqs. (33) and (34)
	if (abs(w) > 0.05) {
		xi0 = sin(w) / w;
	}
	else {
		xi0 = 1.0 - w_squared * (1.0 - w_squared * (1 - w_squared / 42.0) / 20.0) / 6.0;
	}

	
	C_double f0, f1, f2;
	if (c1>0.0001) {
		//if c1 is close to zero, 
		double div = 1.0 / (9.0 * u_squared - w_squared);
		C_double e_pow2Iu((2.0 * u_cos * u_cos - 1.0), 2.0 * u_cos * u_sin);
		C_double e_powminIu(u_cos, -u_sin);
		C_double curlyBrackets(8.0 * u_squared * w_cos, 2.0 * u * (3.0 * u_squared + w_squared) * xi0);
		f0 = div * ((u_squared - w_squared) * e_pow2Iu + e_powminIu * curlyBrackets);
		curlyBrackets.R() = 2.0 * u * w_cos;
		curlyBrackets.I() = (w_squared - 3.0 * u_squared) * xi0;
		f1 = div * (2.0 * u * e_pow2Iu - e_powminIu * curlyBrackets);
		curlyBrackets.R() = w_cos;
		curlyBrackets.I() = 3.0 * u * xi0;
		f2 = div * (e_pow2Iu - e_powminIu * curlyBrackets);

		if (c0 < 0.0) {//eq (34)
			f0 = f0.dagger();
			f1 = -1.0 * f1.dagger();
			f2 = f2.dagger();
		}
	}
	else {
		double c0_squared = c0*c0;
		f0 = C_double(1.0 - c0_squared / 720.0, -c0 * (1.0 - c1 * (1 - c1 / 42.0) / 20.0) / 6.0);
		f1 = C_double(c0 * (1.0 - c1 * (1.0 - 3.0 * c1 / 112.0) / 15.0) / 24.0, 1.0 - c1 * (1.0 - c1 * (1 - c1 / 42.0) / 20.0) / 6.0 - c0_squared / 5040.0);
		f2 = C_double((-1.0 + c1 * (1.0 - c1 * (1.0 - c1 / 56.0) / 30.0) / 12.0 + c0_squared / 20160.0)/2.0, (c0 * (1.0 - c1 * (1.0 - c1 / 48.0) / 21.0) / 60.0)/2.0);
	}

	Q = f1 * Q + f2 * Q_squared;
	Q[0] += f0;
	Q[4] += f0;
	Q[8] += f0;
	return Q;
}

bool IsHermTrLess(su3_mat& mat,bool silence){
	C_double tolerance(0.001,0.001);
	C_double trace = mat.Tr();
	su3_mat matSum(mat - mat.dagger());
	bool out = true;
	if (trace > tolerance) {
		if (!silence) {
			std::cout << "Not hermitian due to trace: " << mat.Tr() << "\n";
		}
		out = false;
	}
	for (int i = 0; i < 9; i++) {
		if (matSum[i] > tolerance) {
			if (!silence) {
				std::cout << "Not hermitian due to entry" << i << "of mat-mat^dagger being" << matSum[i] << "\n";
			}
			out = false;
		}
	}
	return out;
}

std::ostream& operator << (std::ostream& stream, const su3_mat& a) {
	stream << a[0] << " " << a[1] << " " << a[2] << "\n" << a[3] << " " << a[4] << " " << a[5] << "\n" << a[6] << " " << a[7] << " " << a[8];
	return stream;
}
//DDDDDDDDDDDDEEEEEEEEEEEEEEEEELLLLLLLLLLLLLLLLLLEEEEEEEEEEEEEEEETTTTTTTTTTTTTTEEEEEEEEEEEEEEEEEMMMMMMMMMMMMEEEEEEEEEEEEEE
//su3_mat& su3_mat::at()
//{
//	for (int a = 0; a < 3; ++a) {
//		for (int b = a + 1; b < 3; ++b) {
//			double re = mat[(3 * a + b)].R() - mat[(3 * b + a)].R();
//			double im = mat[(3 * a + b)].I() + mat[ (3 * b + a)].I();
//
//			mat[(3* a + b)].R() = 0.5 * re;
//			mat[(3 * a + b)].I() = 0.5 * im;
//
//			mat[(3 * b + a)].R() = -0.5 * re;
//			mat[(3 * b + a)].I() = 0.5 * im;
//		}
//	}
//	double tr = 0.0;
//	for (int cc = 0; cc < 3; ++cc) {
//		tr += mat[(3 * cc + cc)].I();
//	}
//	tr = tr / 3;
//	for (int cc = 0; cc < 3; ++cc) {
//		mat[(3 * cc + cc)].R() = 0.0;
//		mat[(3 * cc + cc)].I() -= tr;
//	}
//	return *this;
//}