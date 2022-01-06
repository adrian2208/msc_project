#pragma once
#include "C_double.h"
#include <iostream>
class su3_mat {
public:
	su3_mat();
	su3_mat(const su3_mat& matrix);

	C_double& operator()(int row, int col);
	C_double& operator[](int i) const;

	su3_mat dagger() const;
	C_double det() const;
	void setToIdentity();

	C_double* mat;
private:
	
	void allocate_matrix();

	

};

su3_mat operator*(const su3_mat& a, const su3_mat& b);
su3_mat operator*(const su3_mat& a, const int b);
su3_mat operator*(const int b, const su3_mat& a);
std::ostream& operator << (std::ostream& stream, const su3_mat& a);