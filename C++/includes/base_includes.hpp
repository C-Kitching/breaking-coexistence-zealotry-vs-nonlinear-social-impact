#pragma once

#ifndef BASE_INCLUDES
#define BASE_INCLUDES

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <utility>
#include <chrono>
#include <omp.h>
#include <filesystem>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <numeric>
#include <Eigen/Dense>

namespace fs = std::filesystem;

// boost multiprecision
// #include <boost/multiprecision/cpp_dec_float.hpp>
// #include <boost/multiprecision/gmp.hpp>
// namespace mp = boost::multiprecision;

// random number generation
#include "../lib/PRNG/pcg_random.hpp"
#include "../lib/PRNG/pcg_extras.hpp"
#include "../lib/PRNG/pcg_uint128.hpp"

typedef std::chrono::high_resolution_clock Clock;
typedef std::mersenne_twister_engine<uint_fast32_t,
	32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253>
	mt19937;

#endif