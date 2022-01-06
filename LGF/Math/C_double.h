#pragma once
#include <iostream>
class C_double {
public:
	C_double();
	C_double(double Real, double Imaginary);
	double& R();
	double& I();
	const double& R() const;
	const double& I() const;

	C_double dagger() const;
	
	
	double Re, Im;
	//C_double operator * (C_double a);
	//C_double operator + (C_double a);
private:


};

C_double operator+(const C_double& a, const C_double& b);
C_double operator-(const C_double& a, const C_double& b);
C_double operator*(const C_double& a, const C_double& b);
C_double operator*(const C_double& a, const int b);
C_double operator*(const int b, const C_double& a);
std::ostream& operator << (std::ostream& stream, const C_double& a);