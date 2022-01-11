#include "C_double.h"

C_double::C_double() {
	Re = 0.0;
	Im = 0.0;
}
C_double::C_double(const double Real = 0.0, const double Imaginary = 0.0) {
	Re = Real;
	Im = Imaginary;
}

double& C_double::R(){
	return Re;
}

double& C_double::I(){
	return Im;
}

const double& C_double::R() const{
	return Re;
}

const double& C_double::I() const{
	return Im;
}

void C_double::operator+=(const C_double& a){
	Re += a.R();
	Im += a.I();
}

void C_double::operator+=(const double a){
	Re += a;
}

void C_double::operator-=(const C_double& a){
	Re -= a.R();
	Im -= a.I();
}

void C_double::operator-=(const double a){
	Re -= a;
}

C_double C_double::dagger() const{
	return C_double(Re,-Im);
}
C_double C_double::timesI() const {
	return C_double(-Im, Re);
}
C_double C_double::timesMinusI() const{
	return C_double(Im,-Re);
}
C_double operator+(const C_double& a, const C_double& b) {
	return C_double(a.R() + b.R(), a.I() + b.I());
}
C_double operator-(const C_double& a, const C_double& b) {
	return C_double(a.R() - b.R(), a.I() - b.I());
}
C_double operator*(const C_double& a, const C_double& b) {
	return C_double(a.R() * b.R() - a.I() * b.I(), a.R() * b.I() + a.I() * b.R());
}
C_double operator*(const C_double& a, const int b){
	return C_double(a.R()*b,a.I()*b);
}
C_double operator*(const int b, const C_double& a){
	return a*b;
}
C_double operator*(const C_double& a, const double b){
	return C_double(a.R()*b,a.I()*b);
}
C_double operator*(const double b, const C_double& a){
	return a*b;
}
std::ostream& operator << (std::ostream& stream, const C_double& a) {
	stream << "(" << a.R();
	if (a.I() == 0) {
		stream << ")";
		return stream;
	}
	if (a.I() > 0) {
		stream << "+";
	}
	stream << a.I() << "i)";
	return stream;
}
//C_double C_double::operator * (C_double a){
//	return C_double(a.Re*Re-a.Im*Im,a.Re*Im+a.Im*Re);
//}
//
//C_double C_double::operator + (C_double a){
//	return C_double(a.Re+Re,a.Im+Im);
//}
