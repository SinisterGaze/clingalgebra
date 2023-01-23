#pragma once
#include "defines.h"
#include "complex.h"


struct dimension
{
	unsigned rows;
	unsigned cols;

	bool operator == (const dimension& rhs) const { return (rows == rhs.rows && cols == rhs.cols); }
	bool operator != (const dimension& rhs) const { return !((*this) == rhs); }

	const std::string str() const { return "(" + std::to_string(rows) + ", " + std::to_string(cols) + ")"; }
	friend std::ostream& operator << (std::ostream& o, const dimension& rhs) { o << rhs.str(); return o; }
};

template<class T>
class Matrix
{
public:
	Matrix();
	Matrix(unsigned rows, unsigned cols);
	Matrix(const Matrix& other);
	Matrix(const std::initializer_list<std::initializer_list<T>>& mat);
	~Matrix();

	unsigned rows() const;
	unsigned cols() const;
	dimension dim() const;

	T& operator () (unsigned row, unsigned col);
	T operator () (unsigned row, unsigned col) const;
	Matrix& operator = (const Matrix& other);
	Matrix& operator = (const std::initializer_list<std::initializer_list<T>>& mat);
	void operator += (const Matrix& rhs);
	void operator -= (const Matrix& rhs);
	Matrix operator + (const Matrix& rhs) const;
	Matrix operator - (const Matrix& rhs) const;
	//Matrix operator * (const Matrix& rhs);
	//Matrix operator * (const T& rhs);
	//void operator *= (const T& rhs);
	//void operator /= (const T& rhs);
	//Matrix operator -();



protected:
	unsigned rows_, cols_;
	T* data_;
};

template <class T>
unsigned Matrix<T>::rows() const { return rows_; }

template <class T>
unsigned Matrix<T>::cols() const { return cols_; }

template <class T>
Matrix<T>::Matrix() {}

template <class T>
Matrix<T>::Matrix(unsigned rows, unsigned cols)
	: rows_(rows), cols_(cols)
{
	if (rows == 0 || cols == 0)
		throw std::invalid_argument("Matrix constructor has 0 size");
	data_ = new T[rows * cols];
}

template <class T>
Matrix<T>::Matrix(const Matrix& other)
	: rows_(other.rows_), cols_(other.cols_)
{
	data_ = new T[rows_ * cols_];
	std::copy(other.data_, other.data_ + rows_ * cols_, data_);
}

template <class T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& mat) : rows_((unsigned)mat.size()), cols_(0)
{
	for (const auto& row : mat)
	{
		unsigned c = (unsigned)row.size();
		if (c > cols_)
			cols_ = c;
	}
	data_ = new T[rows_ * cols_]{ 0 };
	auto it = mat.begin();
	for (unsigned i = 0; i < rows_; ++i, ++it)
		std::copy(it->begin(), it->end(), data_ + cols_ * i);
}

template <class T>
Matrix<T>::~Matrix() { delete data_; }

template <class T>
dimension Matrix<T>::dim() const { return { rows_, cols_ }; }

template <class T>
T& Matrix<T>::operator () (unsigned row, unsigned col)
{
	if (row >= rows_ || col >= cols_)
		throw std::out_of_range("Matrix subscript out of bounds");
	return data_[cols_ * row + col];
}

template <class T>
T Matrix<T>::operator () (unsigned row, unsigned col) const
{
	if (row >= rows_ || col >= cols_)
		throw std::out_of_range("Matrix subscript out of bounds");
	return data_[cols_ * row + col];
}

template <class T>
Matrix<T>& Matrix<T>::operator = (const Matrix& other)
{
	Matrix result(other);
	return result;
}



template<class T>
Matrix<T>& Matrix<T>::operator=(const std::initializer_list<std::initializer_list<T>>& mat)
{
	rows_ = (unsigned)mat.size();
	for (const auto& row : mat)
	{
		unsigned c = (unsigned)row.size();
		if (c > cols_)
			cols_ = c;
	}
	data_ = new T[rows_ * cols_]{ 0 };
	auto it = mat.begin();
	for (unsigned i = 0; i < rows_; ++i, ++it)
		std::copy(it->begin(), it->end(), data_ + cols_ * i);
	return *this;
}

template <class T>
void Matrix<T>::operator += (const Matrix& rhs)
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions does not match for addition");

	for (unsigned i = 0; i < rows_; ++i)
		for (unsigned j = 0; j < cols_; ++j)
			(*this)(i, j) += rhs(i, j);
}

template <class T>
void Matrix<T>::operator -= (const Matrix& rhs)
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions does not match for addition");

	for (unsigned i = 0; i < rows_; ++i)
		for (unsigned j = 0; j < cols_; ++j)
			(*this)(i, j) -= rhs(i, j);
}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix& rhs) const
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions does not match for addition");

	Matrix result(*this);
	result += rhs;
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator - (const Matrix& rhs) const
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions does not match for addition");

	Matrix result(*this);
	result -= rhs;
	return result;
}




typedef Matrix<float> Matrixf;
typedef Matrix<Complexf> Matrixc;


