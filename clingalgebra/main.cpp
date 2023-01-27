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

void test_naive_vs_thread()
{
	for (int rep = 0; rep < 10; rep++)
	{
		int N = 600;
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
		Matrixf D = Matrixf::_thread_mult(A, B);
		end = std::chrono::steady_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time threaded: %lld ms\n", elapsed);

		//print_matrix(C);
		//print_matrix(D);
		cout << (C == D) << endl;
	}
}

void test_thread()
{
	for (int rep = 0; rep < 100; rep++)
	{
		int N = 1000;
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
		Matrixf C = Matrixf::_thread_mult_cached(A, B);
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time threaded mult: %lld ms\n", elapsed);
	}
}

void test_blocking_vs_naive()
{
	for (int rep = 0; rep < 10; rep++)
	{
		int N = 500;
		int R = 800;
		int M = 600;
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
		Matrixf C = Matrixf::_naive_mult(A, B);
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("naive mult: %lld ms\n", elapsed);

		start = std::chrono::steady_clock::now();
		Matrixf D = Matrixf::_fast_mult(A, B);
		end = std::chrono::steady_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("cache blocked mult: %lld ms\n", elapsed);

		cout << (C == D) << "\n";
	}
}

void test_thread_vs_threadcache()
{
	for (int rep = 0; rep < 30; rep++)
	{
		int M = 200;
		int R = 800;
		int N = 500;
		Matrixf A(M, R);
		Matrixf B(R, N);
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < R; j++)
			{
				A(i, j) = (float)(rand() % 50 - 25);
			}
		}
		for (int i = 0; i < R; i++)
		{
			for (int j = 0; j < N; j++)
			{
				B(i, j) = (float)(rand() % 50 - 25);
			}
		}
		printf("N = %d\n", N);
		auto start = std::chrono::steady_clock::now();
		Matrixf C = Matrixf::_thread_mult(A, B);
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time threaded: %lld ms\n", elapsed);

		start = std::chrono::steady_clock::now();
		Matrixf D = Matrixf::_thread_mult_cached(A, B);
		end = std::chrono::steady_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time threaded + cached: %lld ms\n", elapsed);

	}
}

void test_matrix4d_calc()
{

	for (int rep = 0; rep < 100; rep++)
	{
		int N = 4;
		int C = 1'000'000;
		Matrixf T(N, N);
		Matrixf points(N, C);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				T(i, j) = (float)(rand() % 50 - 25);
			}
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < C; j++)
			{
				points(i, j) = (float)(rand() % 50 - 25);
			}
		}

		auto start = std::chrono::steady_clock::now();
		Matrixf newpoints = T * points;
		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("time: %lld ms\n", elapsed);
	}
}

void test_stochastic()
{
	int N = 10;
	Matrixf A(N, N);
	for (int i = 0; i < N; i++)
	{
		Matrixf C(N, N);
		float sum = 0;
		for (int j = 0; j < N; j++)
		{
			float val = (float)(rand() % 50 - 25);
			C(j, i) = val;
			sum += val;
		}
		A += C / sum;
	}
	print_matrix(A);
	for (int i = 0; i < 10; i++)
	{
		A *= A;
		print_matrix(A);
	}
}

int main()
{
	srand((unsigned)time(NULL));
	test_stochastic();


	return 0;
}