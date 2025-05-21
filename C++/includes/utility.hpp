#pragma once

#ifndef UTILITY
#define UTILITY

#include "base_includes.hpp" 


/**
 * @brief Generate equally spaced numbers between 2 bounds with a given amount
 * of required numbers
 * 
 * @tparam T 
 * @param start_in lower bound
 * @param end_in upper bound
 * @param num_in number of numbers to generate
 * @return std::vector<T> 
 */
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if(num == 0) return linspaced;
    if(num == 1){
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i){
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); 

    return linspaced;
}

/**
 * @brief Generate vector of evenly spaced integers from [start, end] inclusive
 * 
 * @tparam T 
 * @param start 
 * @param end 
 * @param N number of points
 * @return std::vector<T> 
 */
template <typename T>
std::vector<T> integer_linspace(T start, T end, std::size_t N)
{
    std::vector<T> result;
    if (N == 0) return result; // No values to generate
    if (N == 1) { // Single value case
        result.push_back(start);
        return result;
    }
    result.reserve(N);
    T step = static_cast<T>(end - start) / (N - 1);
    for (std::size_t i = 0; i < N; ++i) {
        result.push_back(static_cast<T>(std::round(start + i * step)));
    }
    return result;
}

/**
 * @brief Generate all integers [start, end]
 * 
 * @param start 
 * @param end 
 * @return std::vector<int> 
 */
std::vector<int> arange(int start, int end) 
{
    std::vector<int> result;
    if (start <= end) {
        result.reserve(end - start + 1); // Preallocate memory for efficiency
        for (int i = start; i <= end; ++i) {
            result.push_back(i);
        }
    } else {
        result.reserve(start - end + 1); // Preallocate memory for reverse range
        for (int i = start; i >= end; --i) {
            result.push_back(i);
        }
    }
    return result;
}


#endif