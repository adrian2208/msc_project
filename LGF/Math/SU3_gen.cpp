#include "SU3_gen.h"
#include <cmath>
SU3_gen::SU3_gen(){
	su3_mat GM;
	GM.setToZeros();
	C_double one(0.5, 0.0);
	C_double im(0.0, 0.5);
	//1
	GM[1] = GM[3] = one;
	T.push_back(GM);
	//2
	GM[1] = -1.0*im;
	GM[3] = im;
	T.push_back(GM);
	//3
	GM.setToZeros();
	GM[0] = one;
	GM[4] = -1.0 * one;
	T.push_back(GM);
	//4
	GM.setToZeros();
	GM[2] = GM[6] = one;
	T.push_back(GM);
	//5
	GM[2] = -1.0 * im;
	GM[6] = im;
	T.push_back(GM);
	//6
	GM.setToZeros();
	GM[5] = GM[7] = one;
	T.push_back(GM);
	//7
	GM[5] = -1.0 * im;
	GM[7] = im;
	T.push_back(GM);
	//8
	GM.setToZeros();
	GM[0] = GM[4] = 1 / sqrt(3.0) * one;
	GM[8] = -2.0 / sqrt(3.0) * one;
	T.push_back(GM);
}

su3_mat& SU3_gen::operator() (int i) {
	return T[i];
}