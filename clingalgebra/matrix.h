#pragma once
#include "defines.h"
#include "complex.h"

// Helper struct for matrices
struct dimension
{
	unsigned rows;
	unsigned cols;

	bool operator == (const dimension& rhs) const { return (rows == rhs.rows && cols == rhs.cols); }
	bool operator != (const dimension& rhs) const { return !((*this) == rhs); }

	const std::string str() const { return "(" + std::to_string(rows) + ", " + std::to_string(cols) + ")"; }
	friend std::ostream& operator << (std::ostream& o, const dimension& rhs) { o << rhs.str(); return o; }
};

// Template matrix class with unrestricted dimensions
// The template class T is assumed to support the following operations:
// +, -, *, /, +=, -=
// Elements are accessed by m(i, j) = element in i:th row, j:th column
// The copy operator defaults to making deep copies

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
	Matrix operator * (const Matrix& rhs);
	//Matrix operator * (const T& rhs);
	//Matrix operator / (const T& rhs);
	//void operator *= (const T& rhs);
	//void operator /= (const T& rhs);
	//Matrix operator -();
	bool operator == (const Matrix& rhs) const;

	static Matrix _naive_mult(const Matrix& A, const Matrix& B);
	static Matrix _fast_mult(const Matrix& A, const Matrix& B);
	static Matrix _thread_mult(const Matrix& A, const Matrix& B);
	static Matrix _thread_mult_cached(const Matrix& A, const Matrix& B);

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
	cols_ = other.cols_;
	rows_ = other.rows_;
	data_ = new T[rows_ * cols_];
	std::copy(other.data_, other.data_ + rows_ * cols_, data_);
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const std::initializer_list<std::initializer_list<T>>& mat)
{
	rows_ = (unsigned)mat.size();
	cols_ = 0;
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
		throw std::out_of_range("Matrix dimensions incompatible for addition");

	for (unsigned i = 0; i < rows_; ++i)
		for (unsigned j = 0; j < cols_; ++j)
			(*this)(i, j) += rhs(i, j);
}

template <class T>
void Matrix<T>::operator -= (const Matrix& rhs)
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions incompatible for addition");

	for (unsigned i = 0; i < rows_; ++i)
		for (unsigned j = 0; j < cols_; ++j)
			(*this)(i, j) -= rhs(i, j);
}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix& rhs) const
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions incompatible for addition");

	Matrix result(*this);
	result += rhs;
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator - (const Matrix& rhs) const
{
	if (dim() != rhs.dim())
		throw std::out_of_range("Matrix dimensions incompatible for addition");

	Matrix result(*this);
	result -= rhs;
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator * (const Matrix& rhs)
{
	if (cols_ != rhs.rows_)
		throw std::out_of_range("Matrix dimensions incompatible for multiplication");

	unsigned r = rows_;
	unsigned c = rhs.cols_;

	// threading seems to overtake in performance when n ~ 8 * threads
	if (r * c < 8 * 8 * NTHREADS * NTHREADS)
		return _naive_mult(*this, rhs);

	// If matrix is sufficiently large, use threading
	else
	{
		if (r * c < NTHREADS * CACHE_BLOCK * CACHE_BLOCK)
			return _thread_mult(*this, rhs);
		else
			return _thread_mult_cached(*this, rhs);
	}

}

template<class T>
inline bool Matrix<T>::operator==(const Matrix& rhs) const
{
	if (dim() != rhs.dim())
		return false;

	for (unsigned i = 0; i < rows_ * cols_; ++i)
		if (data_[i] != rhs.data_[i])
			return false;

	return true;
}

template<class T>
inline Matrix<T> Matrix<T>::_naive_mult(const Matrix& A, const Matrix& B)
{
	if (A.cols_ != B.rows_)
		throw std::out_of_range("Matrix dimensions incompatible for multiplication");

	unsigned r = A.rows_;
	unsigned c = B.cols_;

	Matrix result(r, c);
	std::fill(result.data_, result.data_ + r * c, (T)0);

	for (unsigned i = 0; i < r; ++i)
		for (unsigned k = 0; k < A.cols_; ++k)
			for (unsigned j = 0; j < c; ++j)
				result.data_[i * c + j] += A.data_[A.cols_ * i + k] * B.data_[B.cols_ * k + j]; //result(i, j) += A(i, k) * B(k, j);

	return result;
}

template<class T>
inline Matrix<T> Matrix<T>::_fast_mult(const Matrix& A, const Matrix& B)
{
	if (A.cols_ != B.rows_)
		throw std::out_of_range("Matrix dimensions incompatible for multiplication");

	unsigned r = A.rows_;
	unsigned c = B.cols_;

	Matrix result(r, c);
	std::fill(result.data_, result.data_ + r * c, (T)0);

	for (int i = 0; i < r; i += CACHE_BLOCK)
		for (int j = 0; j < c; j += CACHE_BLOCK)
			for (int k = 0; k < A.cols_; k += CACHE_BLOCK)
			{
				T* res = &result.data_[i * c + j];
				T* Amul = &A.data_[i * A.cols_ + k];
				int imax = (CACHE_BLOCK < r - i) ? CACHE_BLOCK : r - i;
				for (int i2 = 0; i2 < imax; i2++, res += c, Amul += A.cols_)
				{
					T* Bmul = &B.data_[k * B.cols_ + j];
					int kmax = (CACHE_BLOCK < A.cols_ - k) ? CACHE_BLOCK : A.cols_ - k;
					for (int k2 = 0; k2 < kmax; k2++, Bmul += B.cols_)
					{
						int jmax = (CACHE_BLOCK < c - j) ? CACHE_BLOCK : c - j;
						for (int j2 = 0; j2 < jmax; j2++)
							res[j2] += Amul[k2] * Bmul[j2];
					}
				}
			}

	return result;
}

template<class T>
inline Matrix<T> Matrix<T>::_thread_mult(const Matrix& A, const Matrix& B)
{
	if (A.cols_ != B.rows_)
		throw std::out_of_range("Matrix dimensions incompatible for multiplication");

	unsigned r = A.rows_;
	unsigned c = B.cols_;

	Matrix result(r, c);
	std::fill(result.data_, result.data_ + r * c, (T)0);

	unsigned rowdiv = (unsigned)std::round(std::sqrt(NTHREADS * (float)r / (float)c));
	rowdiv = clingalg::clamp((int)rowdiv, 1, NTHREADS);
	unsigned coldiv = NTHREADS / rowdiv;
	unsigned rblock = r / rowdiv;
	unsigned cblock = c / coldiv;
	rblock = std::max((int)rblock, 1);
	cblock = std::max((int)cblock, 1);

	auto mult_block = [&](unsigned r_idx, unsigned c_idx)
	{
		for (unsigned i = r_idx; i < std::min(r_idx + rblock, r); ++i)
			for (unsigned k = 0; k < A.cols_; k++)
				for (unsigned j = c_idx; j < std::min(c_idx + cblock, c); ++j)
					result.data_[i * c + j] += A.data_[A.cols_ * i + k] * B.data_[B.cols_ * k + j];
	};
	std::vector<std::thread> workers;
	for (unsigned i = 0; i < r; i += rblock)
	{
		for (unsigned j = 0; j < c; j += cblock)
		{
			workers.push_back(std::thread(mult_block, i, j));
		}
	}

	for (auto& worker : workers)
		worker.join();

	return result;
}

template<class T>
inline Matrix<T> Matrix<T>::_thread_mult_cached(const Matrix& A, const Matrix& B)
{
	if (A.cols_ != B.rows_)
		throw std::out_of_range("Matrix dimensions incompatible for multiplication");

	unsigned r = A.rows_;
	unsigned c = B.cols_;

	Matrix result(r, c);
	std::fill(result.data_, result.data_ + r * c, (T)0);

	unsigned rowdiv = (unsigned)std::round(std::sqrt(NTHREADS * (float)r / (float)c));
	rowdiv = clingalg::clamp((int)rowdiv, 1, NTHREADS);
	unsigned coldiv = NTHREADS / rowdiv;
	unsigned rblock = r / rowdiv;
	unsigned cblock = c / coldiv;
	rblock = std::max((int)rblock, 1);
	cblock = std::max((int)cblock, 1);

	auto mult_block = [&](unsigned r_idx, unsigned c_idx)
	{
		int imax = std::min(rblock, r - r_idx);
		for (int i = 0; i < imax; i += CACHE_BLOCK)
		{
			int jmax = std::min(cblock, c - c_idx);
			for (int j = 0; j < jmax; j += CACHE_BLOCK)
			{
				for (int k = 0; k < A.cols_; k += CACHE_BLOCK)
				{
					T* res = &result.data_[(i + r_idx) * c + (j + c_idx)];
					T* Amul = &A.data_[(i + r_idx) * A.cols_ + k];
					int i2max = (CACHE_BLOCK < rblock - i) ? CACHE_BLOCK : rblock - i;
					i2max = std::min(i2max, (int)r - (int)r_idx - (int)i);
					for (unsigned i2 = 0; i2 < i2max; i2++, res += c, Amul += A.cols_)
					{
						T* Bmul = &B.data_[k * B.cols_ + (j + c_idx)];
						int k2max = (CACHE_BLOCK < A.cols_ - k) ? CACHE_BLOCK : A.cols_ - k;
						for (unsigned k2 = 0; k2 < k2max; k2++, Bmul += B.cols_)
						{
							int j2max = (CACHE_BLOCK < cblock - j) ? CACHE_BLOCK : cblock - j;
							j2max = std::min(j2max, (int)c - (int)c_idx - (int)j);
							for (unsigned j2 = 0; j2 < j2max; j2++)
								res[j2] += Amul[k2] * Bmul[j2];
						}
					}
				}
			}
		}
	};
	std::vector<std::thread> workers;
	for (unsigned i = 0; i < r; i += rblock)
	{
		for (unsigned j = 0; j < c; j += cblock)
		{
			workers.push_back(std::thread(mult_block, i, j));
		}
	}
	for (auto& worker : workers)
		worker.join();

	return result;
}






typedef Matrix<float> Matrixf;
typedef Matrix<Complexf> Matrixc;


