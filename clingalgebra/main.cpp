#include "defines.h"
#include "complex.h"
#include "matrix.h"
#include "vector.h"

// rand
#include <cstdlib>

// time measuring
#include <chrono>
#include <time.h>

//printf
#include <cstdio>

// accumulate
#include <numeric>


void print_complex(const Complexf& c)
{
	std::cout << c << std::endl;
}

template<class T>
void print_matrix(const Matrix<T>& m)
{
	unsigned r = m.rows(), c = m.cols();
	for (unsigned i = 0; i < r; ++i)
	{
		std::cout << "[";
		for (unsigned j = 0; j < c; ++j)
		{
			std::cout << m(i, j);
			if (j != c - 1) std::cout << " ";
		}
		std::cout << "]";
		std::cout << "\n";
	}
	std::cout << "\n";
}

using namespace std;

void test_time()
{
	srand(time(NULL));
	int N = 90 * 4;
	const int nreps = 10;
	int times[nreps];
	for (int rep = 0; rep < nreps; rep++)
	{
		Matrixf A(N, N);
		Matrixf B(N, N);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A(i, j) = (float)(rand() % 50 - 25);
				B(i, j) = (float)(rand() % 50 - 25);
			}
		}
		auto start = std::chrono::steady_clock::now();
		Matrixf C = A * B;
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time taken iter %d: %lld ms\n", (int)rep, elapsed);
		times[rep] = (int)elapsed;
	}
	float avg = (float)std::accumulate(times, times + nreps, 0) / (float)nreps;
	printf("average time: %.2f ms\n", avg);
}

void test_mult()
{
	srand(time(NULL));
	for (int rep = 0; rep < 10; rep++)
	{
		int N = 700;
		Matrixf A(N, N);
		Matrixf B(N, N);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A(i, j) = (float)(rand() % 50 - 25);
				B(i, j) = (float)(rand() % 50 - 25);
			}
		}
		auto start = std::chrono::steady_clock::now();
		Matrixf C = Matrixf::_naive_mult(A, B);
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time naive: %lld ms\n", elapsed);

		start = std::chrono::steady_clock::now();
		Matrixf D = A * B;
		end = std::chrono::steady_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time threaded: %lld ms\n", elapsed);

		cout << (C == D) << endl;
	}
}

void test_thread()
{
	srand(time(NULL));
	for (int rep = 0; rep < 30; rep++)
	{
		int N = 500;
		int R = 1500;
		int M = 750;
		Matrixf A(N, R);
		Matrixf B(R, M);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < R; j++)
			{
				A(i, j) = (float)(rand() % 50 - 25);

			}
		}
		for (int i = 0; i < R; i++)
		{
			for (int j = 0; j < M; j++)
			{
				B(i, j) = (float)(rand() % 50 - 25);
			}
		}
		auto start = std::chrono::steady_clock::now();
		Matrixf C = Matrixf::_thread_mult(A, B);
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time naive: %lld ms\n", elapsed);
	}
}



int main()
{
	test_mult();
	return 0;
}