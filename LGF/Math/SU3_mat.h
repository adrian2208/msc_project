#pragma once
#include "C_double.h"
class su3_mat {
public:
	su3_mat();
	su3_mat(const su3_mat& matrix);

	C_double& operator()(int row, int col);
	C_double& operator[](int i) const;

private:
	C_double* mat;
	void allocate_matrix();

	

};