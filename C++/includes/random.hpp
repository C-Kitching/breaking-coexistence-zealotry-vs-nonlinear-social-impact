#pragma once

#ifndef RANDOM_NUMBER_GENERATION
#define RANDOM_NUMBER_GENERATION

#include "base_includes.hpp"

/**
 * @brief Generate a random uniform number
 * 
 * @param rng pseudo random number generator
 * @return double
 */
double draw_random_uniform(pcg32& rng)
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(rng);
}

/**
 * @brief Create a random number generator
 * 
 * @return pcg32 
 */
pcg32 get_random_number_generator()
{
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 source(seed_source);
    return source;
}

#endif