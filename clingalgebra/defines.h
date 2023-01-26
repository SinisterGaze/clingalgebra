#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <thread>
#include <vector>
#include <algorithm>

// Used in matrixmultiplication
#define NTHREADS 12
#define CACHE_BLOCK 16

namespace clingalg
{
	inline int clamp(int v, int lo, int hi)
	{
		return std::max(std::min(v, hi), lo);
	}
}
