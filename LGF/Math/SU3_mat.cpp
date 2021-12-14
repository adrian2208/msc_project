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
