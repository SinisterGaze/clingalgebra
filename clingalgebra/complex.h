#pragma once
#include "defines.h"
template<class T>
class Complex
{
public:
	T a, b;
	Complex() : a(0), b(0) {}
	Complex(T a, T b) : a(a), b(b) {}
	Complex(T a) : a(a), b(0) {}

	T real() const { return a; }
	T imag() const { return b; }
	T arg() const { return T(std::atan2(a, b)); }
	T abs() const { return T(std::sqrt(a * a + b * b)); }
	T norm2() const { return T(a * a + b * b); }

	Complex conj() const { return Complex(a, -b); }

	Complex& operator=(const Complex& rhs) = default;
	Complex& operator=(const T& rhs) { a = rhs; b = 0; return *this; }
	Complex operator + (const Complex& rhs) const { return Complex(this->a + rhs.a, this->b + rhs.b); }
	Complex operator - (const Complex& rhs) const { return Complex(this->a - rhs.a, this->b - rhs.b); }
	Complex operator * (const Complex& rhs) const { return Complex(this->a * rhs.a - this->b * rhs.b, this->a * rhs.b + this->b * rhs.a); }
	Complex operator * (const T& rhs) const { return Complex(this->a * rhs, this->b * rhs); }
	Complex operator / (const T& rhs) const { return Complex(this->a / rhs, this->b / rhs); }
	Complex operator / (const Complex& rhs) { return ((*this) * rhs.conj()) / rhs.norm2(); }
	Complex operator + () { return Complex(+a, +b); }
	Complex operator - () { return Complex(-a, -b); }
	void operator += (const Complex& rhs) { this->a += rhs.a; this->b += rhs.b; }
	void operator -= (const Complex& rhs) { this->a -= rhs.a; this->b -= rhs.b; }
	void operator *= (const T& rhs) { this->a *= rhs; this->b *= rhs; }
	void operator /= (const T& rhs) { this->a /= rhs; this->b /= rhs; }

	const std::string str() const
	{
		std::string sym = (b >= 0 ? "+" : " - ");
		return std::to_string(this->a) + sym + std::to_string(std::abs(this->b)) + "i";
	}
	friend std::ostream& operator << (std::ostream& o, const Complex& rhs) { o << rhs.a << "+" << rhs.b << "i"; return o; }
};

template<class T> inline Complex<T> operator * (const T& lhs, const Complex<T>& rhs) {
	return Complex<T>((T)(lhs * rhs.a), (T)(lhs * rhs.b));
}
template<class T> inline Complex<T> operator * (const int& lhs, const Complex<T>& rhs)
{
	return Complex<T>((T)(lhs * rhs.a), (T)(lhs * rhs.b));
}
template<class T> inline Complex<T> operator * (const float& lhs, const Complex<T>& rhs)
{
	return Complex<T>((T)(lhs * rhs.a), (T)(lhs * rhs.b));
}
template<class T> inline Complex<T> operator * (const double& lhs, const Complex<T>& rhs)
{
	return Complex<T>((T)(lhs * rhs.a), (T)(lhs * rhs.b));
}


template<class T> inline Complex<T> operator * (const Complex<T>& lhs, const T& rhs) {
	return Complex<T>((T)(lhs.a * rhs), (T)(lhs.b * rhs));
}
template<class T> inline Complex<T> operator * (const Complex<T>& lhs, const int& rhs) {
	return Complex<T>((T)(lhs.a * rhs), (T)(lhs.b * rhs));
}
template<class T> inline Complex<T> operator * (const Complex<T>& lhs, const float& rhs) {
	return Complex<T>((T)(lhs.a * rhs), (T)(lhs.b * rhs));
}
template<class T> inline Complex<T> operator * (const Complex<T>& lhs, const double& rhs) {
	return Complex<T>((T)(lhs.a * rhs), (T)(lhs.b * rhs));
}
template<class T> inline Complex<T> operator / (const Complex<T>& lhs, const T& rhs) {
	return Complex<T>((T)(lhs.a / rhs), (T)(lhs.b / rhs));
}
template<class T> inline Complex<T> operator / (const Complex<T>& lhs, const int& rhs) {
	return Complex<T>((T)(lhs.a / rhs), (T)(lhs.b / rhs));
}
template<class T> inline Complex<T> operator / (const Complex<T>& lhs, const float& rhs) {
	return Complex<T>((T)(lhs.a / rhs), (T)(lhs.b / rhs));
}
template<class T> inline Complex<T> operator / (const Complex<T>& lhs, const double& rhs) {
	return Complex<T>((T)(lhs.a / rhs), (T)(lhs.b / rhs));
}

typedef Complex<float> Complexf;