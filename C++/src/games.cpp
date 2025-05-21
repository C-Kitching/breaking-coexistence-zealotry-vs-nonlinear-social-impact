#include "../includes/base_includes.hpp"
#include "../includes/utility.hpp"
#include "../includes/random.hpp"
#include "../includes/progress_bar.hpp"
#include "../includes/file_reading_and_writing.hpp"
#include "../includes/quasi_stationary_dist.hpp"

/**
 * @brief Calculate pi+ for games with absorbing zealots in finite populations
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * @return double pi+
 */
double exp_payoff_plus(const int& N, const int& n, const int& Z,
    const double& alpha)
{
    int S = N-Z;
    double pi_plus = alpha*static_cast<double>(n+Z)/(N-1.) + (1.-alpha)*(S-n)/(N-1.);
    return pi_plus;
}

/**
 * @brief Calculate pi- for games with absorbing zealots in finite populations
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * @return double pi-
 */
double exp_payoff_minus(const int& N, const int& n, const int& Z,
    const double& alpha)
{
    int S = N-Z;
    double pi_minus = (1.-alpha)*static_cast<double>(n+Z)/(N-1.) + alpha*(S-n)/(N-1.);
    return pi_minus;
}

/**
 * @brief Calculate the average payoff for games with absorbing zealots in 
 * finite populations
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * @return double phi
 */
double calculate_avg_payoff(const int& N, const int& n, const int& Z,
    const double& alpha)
{
    double pi_plus = exp_payoff_plus(N, n, Z, alpha);
    double pi_minus = exp_payoff_minus(N, n, Z, alpha);
    int S = N-Z;
    double phi = static_cast<double>(n+Z)/N * pi_plus + static_cast<double>(S-n)/N * pi_minus;
    return phi;
}


/**
 * @brief Calculate T+ analytically for evolutionary games with
 * absorbing zealots
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * 
 * @return double T+
 */
double calculate_T_plus_games(
    const int& N, const int& n, const int& Z, const double& alpha)
{
    int S = N-Z;
    double pi_plus = exp_payoff_plus(N, n, Z, alpha);
    // double phi = calculate_avg_payoff(N, n, Z, alpha);
    // double Tplus = (S-n)*static_cast<double>(n+Z)/(N-1.)*pi_plus/phi;
    double Tplus = (S-n)*static_cast<double>(n+Z)/(N-1.)*pi_plus;
    return Tplus; 
}

/**
 * @brief Calculate T- analytically for evolutionary games with
 * absorbing zealots
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * 
 * @return double T-
 */
double calculate_T_minus_games(
    const int& N, const int& n, const int& Z, const double& alpha)
{
    int S = N-Z;
    double pi_minus = exp_payoff_minus(N, n, Z, alpha);
    // double phi = calculate_avg_payoff(N, n, Z, alpha);
    // double Tminus = n*static_cast<double>(S-n)/(N-1.) * pi_minus/phi;
    double Tminus = n*static_cast<double>(S-n)/(N-1.) * pi_minus;
    return Tminus; 
}

/**
 * @brief Calculate gamma_k for all k=[0, N] for games with absorbing zealots
 * in finite populations
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_gammas_games(
    const int& N, const int& Z, const double& alpha)
{
    std::vector<double> gamma_ks(N+1);
    double Tplus, Tminus;

    for (int n{0}; n <= N; n++){
        Tplus = calculate_T_plus_games(N, n, Z, alpha);
        Tminus = calculate_T_minus_games(N, n, Z, alpha);
        gamma_ks[n] = Tminus/Tplus;
    }

    return gamma_ks;
}

/**
 * @brief Calculate T_n+ analytically for all n in range [0, N] for games
 * with zealots in finite populations
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_all_Tplus_games(
    const int& N, const int& Z, const double& alpha)
{
    std::vector<double> res(N);

    for (int n{0}; n <= N; n++)
        res[n] = calculate_T_plus_games(N, n, Z, alpha);

    return res;
}

/**
 * @brief Calculate T_n- analytically for all n in range [0, N] for games
 * with zealots in finite populations
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param alpha payoff param
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_all_Tminus_games(
    const int& N, const int& Z, const double& alpha)
{
    std::vector<double> res(N);

    for (int n{0}; n <= N; n++)
        res[n] = calculate_T_minus_games(N, n, Z, alpha);

    return res;
}

/**
 * @brief Calculate the fixation time analytically for different z values
 * for games, then write the results to file
 * 
 */
void fixation_time_analytic_against_z_games(const int& N, const double& alpha)
{
    std::stringstream alpha_ss;
    alpha_ss << "alpha=" << alpha;
    int n = 0;
    std::vector<int> Z_store = arange(1, N-1);

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "games", "analytical_fixation_times_against_z", alpha_ss.str()};

    // results containers
    std::vector<double> fixation_times(Z_store.size());

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = Z_store.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // loop over all z
    for (int i{0}; i < static_cast<int>(Z_store.size()); i++){

        double Z = Z_store[i];
        std::vector<double> gammas = calculate_gammas_games(N, Z, alpha);
        std::vector<double> Tplus = calculate_all_Tplus_games(N, Z, alpha);

        double sum1{0}, sum2{0};
        for (int k{n}; k <= N-Z-1; k++){

            double prod1{1};
            for (int m{1}; m <= k; m++) prod1 *= gammas[m];
            sum1 += prod1;

            double inner_sum2{0};
            for (int l{1}; l <= k; l++){

                double prod2{1};
                for (int m{l+1}; m<= k; m++) prod2 *= gammas[m];
                inner_sum2 += prod2/Tplus[l];

            }
            sum2 += inner_sum2;

        }

        fixation_times[i] = sum1/Tplus[0] + sum2;
    
        progress++;
    }
    progress_thread.join();

    // prepare results for writing
    std::vector<std::vector<double>> res(Z_store.size(), std::vector<double>(2));
    for (size_t i = 0; i < Z_store.size(); i++)
        res[i] = {static_cast<double>(Z_store[i])/N, fixation_times[i]};

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "analytical_fixation_times");

}

/**
 * @brief Calculate the average fixation time in simulation over multiple
 * trajectories in the games model with finite populations for different 
 * values of z and then write the data to file
 * 
 */
void average_fixation_times_against_z_games(const int& N, const double& alpha)
{
    std::stringstream alpha_ss;
    alpha_ss << "alpha=" << alpha;
    int n = 0;
    std::vector<int> Z_store = integer_linspace(1, N-1, 25);
    double T = 1e5;
    int number_of_sims = 1e3;

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "games", "average_fixation_times_against_z", alpha_ss.str()};

    // store results
    std::vector<std::atomic<double>> average_fixation_times(Z_store.size());
    std::vector<int> count_fixation_reached(Z_store.size());

    // lambda function to update atomic double
    auto update_atomic_double = [](std::atomic<double>& atom, double val){
        double old_val = atom.load(std::memory_order_relaxed);
        double new_val;
        do {
            new_val = old_val + val;
        } while (!atom.compare_exchange_weak(old_val, new_val, 
                                            std::memory_order_relaxed,
                                            std::memory_order_relaxed));
    };

    // random number generator
    pcg32 source = get_random_number_generator();

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = Z_store.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // loop over all z
    for (size_t i{0}; i < Z_store.size(); i++){
        double Z = Z_store[i];
        int S = N-Z;

        // track the number of sims that reach fixation
        std::atomic<int> count_fixation{0};

        // split over cores
        #pragma omp parallel
        {
            // initialise thread local PRNG and jump ahead
            thread_local pcg32 rng = source;
            rng.advance(omp_get_thread_num() * 10000);

            // loop over all sims
            #pragma omp for
            for (int j = 0; j < number_of_sims; j++){

                // initialise system
                double num_up_nodes = n;
                double Tplus = calculate_T_plus_games(N, num_up_nodes, Z, alpha);
                double Tminus = calculate_T_minus_games(N, num_up_nodes, Z, alpha);
                double total_rate = Tplus + Tminus;

                // run for time T
                double t{0};
                double r1, delta;
                while (t < T){

                    // determine update time
                    r1 = draw_random_uniform(rng);
                    delta = -1.*log(r1)/total_rate;
                    t += delta;

                    // flip node
                    if (draw_random_uniform(rng) <= Tplus/total_rate) num_up_nodes += 1;
                    else num_up_nodes -= 1;

                    // recalculate rates
                    Tplus = calculate_T_plus_games(N, num_up_nodes, Z, alpha);
                    Tminus = calculate_T_minus_games(N, num_up_nodes, Z, alpha);
                    total_rate = Tplus + Tminus;

                    // absorbing state
                    if (num_up_nodes == S){
                        count_fixation++;
                        update_atomic_double(average_fixation_times[i], t);  // accumulate fixation times
                        break;
                    }
                }
            }
        }

        // store number of sims to reach fixation
        count_fixation_reached[i] = count_fixation.load();

        progress++;
    }
    progress_thread.join();

    // prepare results for writing
    std::vector<std::vector<double>> res(Z_store.size(), std::vector<double>(2));
    for (size_t i = 0; i < Z_store.size(); i++)
        res[i] = {static_cast<double>(Z_store[i])/N, average_fixation_times[i]/count_fixation_reached[i]};

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "average_fixation_times");

}

// /**
//  * @brief For a variety of alpha values determine the Z value at which the 
//  * recursive solution to the quasi stationary distribution is peaked at the
//  * upper boundary
//  * 
//  * @param N size of system
//  */
// void quasi_critical_z_games(const int& N)
// {
//     // set up folders for writing data
//     std::vector<std::string> folder_names = {
//         "games", "quasi_stationary_distribution"};

//     std::vector<double> alpha_values = linspace(0.01, 0.99, 30); 
//     std::vector<double> alpha_store; 
//     std::vector<double> corresponding_Z; 

//     // numerical integration params
//     double delta = 0.001;
//     int max_iter = 1e5;
//     double tolerance = 1e-8;

//     // Initialize atomic progress counter
//     std::atomic<int> progress(0);
//     int total_iterations = alpha_values.size();

//     // Start progress display thread
//     std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

//     // loop over q values
//     for (const auto& alpha : alpha_values) {

//         bool found = false;

//         // loop over Z in [1, N-1] (i.e. always both susceptible and zealot)
//         for (int Z = 1; Z <= N-1; Z++) {

//             // initailsie flat distribution
//             int size = N - Z;
//             Eigen::VectorXd q_recursive = Eigen::VectorXd::Ones(size);
//             q_recursive /= q_recursive.sum();

//             // pre calculate all variables
//             std::vector<double> all_Tplus = calculate_all_Tplus_games(N, Z, alpha);
//             std::vector<double> all_Tminus = calculate_all_Tminus_games(N, Z, alpha);
//             Eigen::MatrixXd A = calculate_A_matrix(all_Tplus, all_Tminus, N, Z);
        
//             // Perform numerical integration
//             explicit_euler(A, q_recursive, all_Tplus[size-1], delta, max_iter, tolerance);

//             // Determine the index of the maximum element.
//             auto max_it = std::max_element(q_recursive.begin(), q_recursive.end());
//             int max_index = std::distance(q_recursive.begin(), max_it);

//             // Check if the maximum is at the upper boundary.
//             if (max_index == size - 1) {
//                 alpha_store.push_back(alpha);
//                 corresponding_Z.push_back(Z);
//                 found = true;
//                 break;  // exit the Z loop once the condition is met.
//             }
//         }

//         // Record -1 if no Z produces upper peal
//         if (!found) {
//             alpha_store.push_back(alpha);
//             corresponding_Z.push_back(-1);  // or some flag value indicating no solution.
//         }

//         progress++;

//     }
//     progress_thread.join();

//     // Write to file
//     std::vector<std::vector<double>> res;
//     for (size_t i{0}; i < alpha_store.size(); i++)
//         res.push_back({alpha_store[i], corresponding_Z[i]});

//     // write data to file
//     write_vectors_to_file_rowwise(res, folder_names, "quasi_critical_Z_varying_alpha");

// }