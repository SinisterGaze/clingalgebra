#pragma once
#include "defines.h"
#include "matrix.h"

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

typedef Vector<float> Vectorf;
typedef Vector<Complexf> Vectorc;


template <class T>
class rVector : public Matrix<T>
{
public:
	rVector(unsigned cols) : Matrix<T>(1, cols) {}

	rVector(const std::initializer_list<T>& vec) : Matrix<T>({ vec }) {}

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

typedef rVector<float> rVectorf;
typedef rVector<Complexf> rVectorc;