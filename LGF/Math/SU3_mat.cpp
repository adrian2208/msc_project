#include "SU3_mat.h"
su3_mat::su3_mat(){
	allocate_matrix();
}
su3_mat::su3_mat(const su3_mat& matrix) {
	allocate_matrix();
	for (int i = 0; i < 9; i++) {
		mat[i] = matrix[i];
	}
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

C_double su3_mat::det() const{	
	return mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]) - mat[1] * (mat[3] * mat[8] - mat[5] * mat[6]) + mat[2] * (mat[3] * mat[7] - mat[4] * mat[6]);
}

void su3_mat::setToIdentity(){
	mat[0].Re = mat[4].Re = mat[8].Re = 1.0;
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

std::ostream& operator << (std::ostream& stream, const su3_mat& a) {
	stream << a[0] << " " << a[1] << " " << a[2] << "\n" << a[3] << " " << a[4] << " " << a[5] << "\n" << a[6] << " " << a[7] << " " << a[8];
	return stream;
}