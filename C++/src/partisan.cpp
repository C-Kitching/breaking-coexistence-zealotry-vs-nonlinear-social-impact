#include "../includes/base_includes.hpp"
#include "../includes/utility.hpp"
#include "../includes/random.hpp"
#include "../includes/progress_bar.hpp"
#include "../includes/file_reading_and_writing.hpp"

/**
 * @brief Calculate the average fixation time in simulations for the partisan
 * voter model in finite populations with zealots
 * 
 * @param epsilon noise rate
 * @param N population size
 * 
 */
void average_fixation_times_against_z_partisan(const double& epsilon,
    const double& N)
{
    // Precompute constants
    const double N_inv = 1./(N-1.);
    const double half = 0.5;

    // File set up
    std::stringstream alpha_ss;
    alpha_ss << "epsilon=" << epsilon;
    std::vector<int> Z_store = integer_linspace(1, static_cast<int>(N/2.), 25);

    // simulation params
    double T = 1e7;
    const int number_of_sims = 1e4;

    // Set up folders for writing data
    std::vector<std::string> folder_names = {
        "partisan", "average_fixation_times_against_z", alpha_ss.str()
    };

    // Store results
    std::vector<double> average_fixation_times(Z_store.size(), 0.0);
    std::vector<int> count_fixation_reached(Z_store.size(), 0);

    // Random number generator
    pcg32 source = get_random_number_generator();

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    const int total_iterations = Z_store.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // Parallelize outer loop over Z_store
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < Z_store.size(); i++) {

        double Z = Z_store[i];
        const int S = N - Z;
        // Equal number of agents preferring + or -
        const int plus_pref = S / 2;
        const int minus_pref = S - plus_pref;

        // Thread-local RNG
        thread_local pcg32 rng = source;
        rng.advance(omp_get_thread_num() * 10000 + i);

        // Thread-local variables
        double fixation_time_sum = 0.0;
        int count_fixation = 0;

        // Loop over all simulations
        for (int j = 0; j < number_of_sims; j++) {

            // Initialize system
            int nPP = 0; // x_+^+ i.e. agents in state + who prefer +
            int nPM = 0; // x_+^- i.e. agents in state + who prefer -
            int nMP = plus_pref;
            int nMM = minus_pref;

            // Run for time T
            double t = 0.0;
            double Tplus_pluspreference, Tminus_pluspreference;
            double Tplus_minuspreference, Tminus_minuspreference;
            double total_rate;

            // Initial rates
            Tplus_pluspreference = nMP * (nPP + nPM + Z) * N_inv * (1 + epsilon) * half;
            Tminus_pluspreference = nPP * (nMP + nMM) * N_inv * (1 - epsilon) * half;
            Tplus_minuspreference = nPM * (nMP + nMM) * N_inv * (1 + epsilon) * half;
            Tminus_minuspreference = nMM * (nPP + nPM + Z) * N_inv * (1 - epsilon) * half;
            total_rate = Tplus_pluspreference + Tminus_pluspreference + Tplus_minuspreference + Tminus_minuspreference;

            while (t < T) {
                // Determine update time
                double r1 = draw_random_uniform(rng);
                double delta = -std::log(r1) / total_rate;
                t += delta;

                // Perform event
                double r2 = draw_random_uniform(rng) * total_rate;
                if (r2 <= Tplus_pluspreference) {
                    nPP += 1;
                    nMP -= 1;
                } else if (r2 <= Tplus_pluspreference + Tminus_pluspreference) {
                    nPP -= 1;
                    nMP += 1;
                } else if (r2 <= Tplus_pluspreference + Tminus_pluspreference + Tplus_minuspreference) {
                    nMM += 1;
                    nPM -= 1;
                } else {
                    nMM -= 1;
                    nPM += 1;
                }

                // Check for absorbing state
                if (nPP + nPM == S) {
                    count_fixation++;
                    fixation_time_sum += t;
                    break;
                }

                // Recalculate rates
                Tplus_pluspreference = nMP * (nPP + nPM + Z) * N_inv * (1 + epsilon) * half;
                Tminus_pluspreference = nPP * (nMP + nMM) * N_inv * (1 - epsilon) * half;
                Tplus_minuspreference = nPM * (nMP + nMM) * N_inv * (1 + epsilon) * half;
                Tminus_minuspreference = nMM * (nPP + nPM + Z) * N_inv * (1 - epsilon) * half;
                total_rate = Tplus_pluspreference + Tminus_pluspreference + Tplus_minuspreference + Tminus_minuspreference;

                // Handle possible zero total_rate to prevent division by zero
                if (total_rate <= 0.0) {
                    break; // Cannot proceed further
                }
            }
        }

        // Store results
        average_fixation_times[i] = (count_fixation > 0) ? (fixation_time_sum / count_fixation) : 0.0;
        count_fixation_reached[i] = count_fixation;

        progress++;
    }
    progress_thread.join();

    // Prepare results for writing
    std::vector<std::vector<double>> res(Z_store.size(), std::vector<double>(2));
    for (size_t i = 0; i < Z_store.size(); i++){
        res[i] = {
            Z_store[i] / N,
            average_fixation_times[i]
        };
    }

    // Write data to file
    write_vectors_to_file_rowwise(res, folder_names, "average_fixation_times");
}