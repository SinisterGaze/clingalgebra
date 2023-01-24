#pragma once
#include "defines.h"
#include "matrix.h"

// Template class for "regular" column vectors
// i.e. (n x 1)-matrices
template <class T>
class Vector : public Matrix<T>
{
public:
	Vector() { this->cols_ = 1; }

	Vector(unsigned rows) : Matrix<T>(rows, 1) {}

	Vector(const std::initializer_list<T>& vec) : Matrix<T>((unsigned)vec.size(), 1)
	{
		std::copy(vec.begin(), vec.end(), this->data_);
	}

	Vector& operator = (const std::initializer_list<T>& vec)
	{
		this->rows_ = (unsigned)vec.size();
		this->cols_ = 1;
		this->data_ = new T[this->rows_ * this->cols_];
		std::copy(vec.begin(), vec.end(), this->data_);
		return *this;
	}

	T& operator () (unsigned row)
	{
		if (row >= Matrix<T>::rows_)
			throw std::out_of_range("Vector subscript out of bounds");
		return Matrix<T>::data_[row];
	}

	T operator () (unsigned row) const
	{
		if (row >= Matrix<T>::rows_)
			throw std::out_of_range("Vector subscript out of bounds");
		return Matrix<T>::data_[row];
	}
};



// Template class for row vectors
// i.e. (1 x n)-matrices
template <class T>
class rVector : public Matrix<T>
{
public:
	rVector(unsigned cols) : Matrix<T>(1, cols) {}

	rVector(const std::initializer_list<T>& vec) : Matrix<T>({ vec }) {}

	rVector& operator = (const std::initializer_list<T>& vec)
	{
		this->rows_ = 1;
		this->cols_ = vec.size();
		this->data_ = new T[this->rows_ * this->cols_];
		std::copy(vec.begin(), vec.end(), this->data_);
		return *this;
	}

	T& operator () (unsigned col)
	{
		if (col >= Matrix<T>::cols_)
			throw std::out_of_range("Vector subscript out of bounds");
		return Matrix<T>::data_[col];
	}

	T operator () (unsigned col) const
	{
		if (col >= Matrix<T>::cols_)
			throw std::out_of_range("Vector subscript out of bounds");
		return Matrix<T>::data_[col];
	}
};

typedef Vector<float> Vectorf;
typedef Vector<Complexf> Vectorc;

typedef rVector<float> rVectorf;
typedef rVector<Complexf> rVectorc;