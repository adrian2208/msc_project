#pragma once
#include "C_double.h"
#include <iostream>
#include <vector>


class su3_mat {
public:
	su3_mat();
	su3_mat(const su3_mat& matrix);
	~su3_mat();

	C_double& operator()(int row, int col);
	C_double& operator[](int i) const;
	const su3_mat& operator=(const su3_mat& a);

	su3_mat dagger() const;
	su3_mat dagger();
	su3_mat timesI() const;
	su3_mat timesMinusI() const;
	C_double det() const;
	C_double Tr() const;
	double ReTr() const;
	su3_mat& at();//!!!!!!!!!!!!!!!!!!!!!!DELETE MEEEEE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	void setToIdentity();
	void setToZeros();
	C_double* getMemPointer();
	
private:
	C_double* mat;
	void allocate_matrix();

	

};

su3_mat operator*(const su3_mat& a, const su3_mat& b);
su3_mat operator+(const su3_mat& a, const su3_mat& b);
su3_mat operator-(const su3_mat& a, const su3_mat& b);

su3_mat operator*(const su3_mat& a, const int b);
su3_mat operator*(const int b, const su3_mat& a);

su3_mat operator*(const su3_mat& a, const double b);
su3_mat operator*(const double b, const su3_mat& a);

su3_mat operator*(const su3_mat& a, const C_double b);
su3_mat operator*(const C_double b,const su3_mat& a);
/// <summary>
/// This function will exponentiate any 3x3 hermitian, traceless matrix
/// according to the procedure outlined in the paper
/// "Analytic smearing of SU(3) link variables in lattice QCD" by Colin Morningstar
/// and Mike Peardon
/// </summary>
/// <param name="Q"></param>
su3_mat HermTrLessExp(su3_mat& iP);

bool IsHermTrLess(su3_mat& mat, bool silence= false);
bool IsAntiHermTrLess(su3_mat& mat, bool silence = false);
bool isSpecialUnitary(su3_mat& mat, bool silence = false);

std::ostream& operator << (std::ostream& stream, const su3_mat& a);