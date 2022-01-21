#pragma once
#include <iostream>
class C_double {
public:
	C_double();
	C_double(double Real, double Imaginary);

	inline double& R() {return Re;}
	inline double& I() {return Im;}
	inline const double& R() const {return Re;}
	inline const double& I() const {return Im;}
	void operator +=(const C_double& a);
	void operator +=(const double a);
	void operator -=(const C_double& a);
	void operator -=(const double a);
	//C_double dagger() const;
	//C_double timesI() const;
	//C_double timesMinusI() const;
	inline C_double dagger() const {
		return C_double(Re, -Im);
	}
	inline C_double timesI() const {
		return C_double(-Im, Re);
	}
	inline C_double timesMinusI() const {
		return C_double(Im, -Re);
	}
	
	double Re, Im;
	//C_double operator * (C_double a);
	//C_double operator + (C_double a);
private:


};

//C_double operator+(const C_double& a, const C_double& b);
//C_double operator-(const C_double& a, const C_double& b);
//C_double operator*(const C_double& a, const C_double& b);
inline C_double operator+(const C_double& a, const C_double& b) {
	return C_double(a.Re + b.R(), a.I() + b.I());
}
inline C_double operator-(const C_double& a, const C_double& b) {
	return C_double(a.R() - b.R(), a.I() - b.I());
}
inline C_double operator*(const C_double& a, const C_double& b) {
	return C_double(a.R() * b.R() - a.I() * b.I(), a.R() * b.I() + a.I() * b.R());
}

C_double operator+(const C_double& a, const double b);
C_double operator+(const double b, const C_double& a);
C_double operator*(const C_double& a, const int b);
C_double operator*(const int b, const C_double& a);
C_double operator*(const C_double& a, const double b);
C_double operator*(const double b, const C_double& a);
bool operator<(const C_double& a, const C_double& b);
bool operator>(const C_double& a, const C_double& b);
bool operator>(const C_double& a, const double b);
bool operator>(const double b, const C_double& a);
bool operator<(const C_double& a, const double b);
bool operator<(const double b, const C_double& a);
std::ostream& operator << (std::ostream& stream, const C_double& a);