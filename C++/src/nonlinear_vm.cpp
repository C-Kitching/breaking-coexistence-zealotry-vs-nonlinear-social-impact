#include "../includes/base_includes.hpp"
#include "../includes/utility.hpp"
#include "../includes/random.hpp"
#include "../includes/progress_bar.hpp"
#include "../includes/file_reading_and_writing.hpp"
#include "../includes/quasi_stationary_dist.hpp"

/**
 * @brief Evolve x in time for a single realsiation of the nonlinear voter 
 * model and write the data to a file
 * 
 */
void evolve_x_in_time_single_realisation_nonlinear_vm()
{
    double q = 3;
    std::stringstream qss;
    qss << "q=" << q;
    double z = 0.3;
    std::stringstream zss;
    zss << "z=" << z;
    int N = 100;
    int Z = z*N; 
    int S = N-Z;
    int n = 1;
    double T = 1e3;

    std::vector<std::string> folder_names = {
        "nonlinear_vm", "x_with_time_single_realisation", qss.str(), zss.str()};

    // random number generator
    pcg32 rng = get_random_number_generator();

    double samples = 1000;

    // data arrays
    std::vector<double> x_store(samples+1);
    std::vector<double> time = linspace(0., T, samples+1);

    // initialise system
    double num_up_nodes = n;
    int sampling_counter = 0;
    double Tplus = (S-num_up_nodes)*std::pow((num_up_nodes+Z)/(N-1), q);
    double Tminus = num_up_nodes*std::pow((S-num_up_nodes)/(N-1), q);
    double total_rate = Tplus + Tminus;

    // run for time T
    double t{0};
    double r1, delta;
    while (t < T){

        // determine update time
        r1 = draw_random_uniform(rng);
        delta = -1.*log(r1)/total_rate;
        t += delta;

        // store data
        if (t >= time[sampling_counter]){
            
            // back fill previous times
            for (int j{sampling_counter}; j <= samples; j++){
                if (time[j] <= t){
                    x_store[sampling_counter] = num_up_nodes/N;
                    sampling_counter++;
                } else break;
            }
        }

        // flip node
        if (draw_random_uniform(rng) <= Tplus/total_rate) num_up_nodes += 1;
        else num_up_nodes -= 1;

        // recalculate rates
        Tplus = (S-num_up_nodes)*std::pow((num_up_nodes+Z)/(N-1), q);
        Tminus = num_up_nodes*std::pow((S-num_up_nodes)/(N-1), q);
        total_rate = Tplus + Tminus;

        // absorbing state
        if (num_up_nodes == S){

            // fill remaining times
            while (sampling_counter <= samples){
                x_store[sampling_counter] = num_up_nodes/N;
                sampling_counter++;
            }

            break;
        }
    }

    // prepare results for writing
    std::vector<std::vector<double>> res(time.size(), std::vector<double>(2));
    for (size_t i = 0; i < time.size(); i++)
        res[i] = {time[i], x_store[i]};

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "x_with_time");

}

/**
 * @brief Calculate the fixation time in simulation for a single realisation 
 * in the nonlinear voter model for different values of z
 * 
 */
void fixation_times_against_z_single_realisation_nonlinear_vm()
{
    double q = 1.5;
    std::stringstream qss;
    qss << "q=" << q;
    std::vector<double> z_store = linspace(0.01, 0.99, 1000);
    int N = 100;
    int n = 1;
    double T = 1e3;

    // store results
    std::vector<double> fixation_times(z_store.size());

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "fixation_times_against_z_single_realisation", qss.str()};

    // random number generator
    pcg32 rng = get_random_number_generator();

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = z_store.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // loop over all z
    for (size_t i{0}; i < z_store.size(); i++){
        double z = z_store[i];
        int Z = z*N; 
        int S = N-Z;

        // initialise system
        double num_up_nodes = n;
        double Tplus = (S-num_up_nodes)*std::pow((num_up_nodes+Z)/(N-1), q);
        double Tminus = num_up_nodes*std::pow((S-num_up_nodes)/(N-1), q);
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
            Tplus = (S-num_up_nodes)*std::pow((num_up_nodes+Z)/(N-1), q);
            Tminus = num_up_nodes*std::pow((S-num_up_nodes)/(N-1), q);
            total_rate = Tplus + Tminus;

            // absorbing state
            if (num_up_nodes == S) break;
        }

        // store results
        fixation_times[i] = t;

        progress++;
    }
    progress_thread.join();

    // prepare results for writing
    std::vector<std::vector<double>> res(z_store.size(), std::vector<double>(2));
    for (size_t i = 0; i < z_store.size(); i++)
        res[i] = {z_store[i], fixation_times[i]};

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "fixation_times");

}

/**
 * @brief 
 * 
 *
 * absorbing zealots
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param q
 * @return double T+
 */
double calculate_T_plus_nonlinear_vm(
    const int& N, const int& n, const int& Z, const double& q)
{
    int S = N-Z;
    double Tplus = (S-n)*std::pow((n+Z)/(N-1.), q);
    return Tplus; 
}

/**
 * @brief Calculate T- analytically for the nonlinear voter model with
 * absorbing zealots
 * 
 * @param N population size
 * @param n number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param q
 * @return double T-
 */
double calculate_T_minus_nonlinear_vm(
    const int& N, const int& n, const int& Z, const double& q)
{
    int S = N-Z;
    double Tminus = n*std::pow((S-n)/(N-1.), q);
    return Tminus; 
}

/**
 * @brief Calculate gamma_k for all k=[0, N] in the nonlinear voter model
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param q
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_gammas_nonlinear_vm(
    const int& N, const int& Z, const double& q)
{
    std::vector<double> gamma_ks(N+1);
    double Tplus, Tminus;

    for (int n{0}; n <= N-Z; n++){
        Tplus = calculate_T_plus_nonlinear_vm(N, n, Z, q);
        Tminus = calculate_T_minus_nonlinear_vm(N, n, Z, q);
        gamma_ks[n] = Tminus/Tplus;
    }

    return gamma_ks;
}

/**
 * @brief Calculate T_n+ analytically for all n in range [0, N] for the
 * nonlinear voter model
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param q
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_all_Tplus_nonlinear_vm(
    const int& N, const int& Z, const double& q)
{
    std::vector<double> res(N-Z+1);

    for (int n{0}; n <= N-Z; n++)
        res[n] = calculate_T_plus_nonlinear_vm(N, n, Z, q);

    return res;
}

/**
 * @brief Calculate T_n- analytically for all n in range [0, N] for the
 * nonlinear voter model
 * 
 * @param N population size
 * @param Z number of +1 zealots
 * @param q
 * 
 * @return std::vector<double> 
 */
std::vector<double> calculate_all_Tminus_nonlinear_vm(
    const int& N, const int& Z, const double& q)
{
    std::vector<double> res(N-Z+1);

    for (int n{0}; n <= N-Z; n++)
        res[n] = calculate_T_minus_nonlinear_vm(N, n, Z, q);

    return res;
}

/**
 * @brief Calculate the fixation time analytically for different z valus in the
 * nonlinear voter model, then write the results to file
 * 
 */
void fixation_time_analytic_against_z_nonlinear_vm()
{

    double q = 3;
    std::stringstream qss;
    qss << "q=" << q;
    int N = 100;
    int n = 0;
    std::vector<int> Z_store = arange(1, N-1);

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "analytical_fixation_times_against_z", qss.str()};

    // results containers
    std::vector<double> fixation_times(Z_store.size());

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = Z_store.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // loop over all z
    for (int i{0}; i < static_cast<int>(Z_store.size()); i++){

        int Z = Z_store[i];
        std::vector<double> gammas = calculate_gammas_nonlinear_vm(N, Z, q);
        std::vector<double> Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);

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
 * trajectories in the nonlinear voter model for different values of z
 * and then write the data to file
 * 
 */
void average_fixation_times_against_z_nonlinear_vm()
{
    double q = 0.8;
    std::stringstream qss;
    qss << "q=" << q;
    int N = 100;
    int n = 0;
    std::vector<int> Z_store = integer_linspace(1, N-1, 25);
    double T = 1e3;
    int number_of_sims = 1e3;

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "average_fixation_times_against_z", qss.str()};

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
        int Z = Z_store[i];
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
                double Tplus = calculate_T_plus_nonlinear_vm(N, num_up_nodes, Z, q);
                double Tminus = calculate_T_minus_nonlinear_vm(N, num_up_nodes, Z, q);
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
                    Tplus = calculate_T_plus_nonlinear_vm(N, num_up_nodes, Z, q);
                    Tminus = calculate_T_minus_nonlinear_vm(N, num_up_nodes, Z, q);
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

/**
 * @brief Calculate a single fixation time analytically for a specific
 * set of parameters. This is for nonlinear voter model with absorbing 
 * zealots in finite populations
 * 
 * @param N population size
 * @param n initial number of +1 susceptibles
 * @param Z number of +1 zealots
 * @param q 
 * @return double t_0
 */
double calculate_specific_fixation_time_nonlinear_vm(
    const int& N, const int& n, const int& Z, const double& q)
{
    std::vector<double> gammas = calculate_gammas_nonlinear_vm(N, Z, q);
    std::vector<double> Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);

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

    double fixation_time = sum1/Tplus[0] + sum2;

    return fixation_time;
}

/**
 * @brief Calculate the quasi-stationary distribution in simulations for the
 * nonlinear voter model in finite populations with absorbing zealots,
 * using a cloning (resampling) algorithm.
 *  
 * @param q nonlinear vm param
 * @param N system size
 * @param Z number of zealots
 */
void quasi_stationary_distribution_simulation_nonlinear_vm(const double& q, const int& N, const int& Z)
{
    std::stringstream qss;
    qss << "q=" << q;
    std::stringstream Zss;
    Zss << "Z=" << Z;
    int S = N - Z; // absorbing state: S (number of susceptible nodes)
    int number_of_sims = 1e4;

    // time to run the simulations needs to be tuned for the specific situation
    double fixation_time = calculate_specific_fixation_time_nonlinear_vm(N, 0, Z, q);
    double T = 2 * fixation_time; // q<1
    // double T = 1000;  // q>1

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "quasi_stationary_distribution", qss.str(), Zss.str()
    };

    // ---------------------------------------
    // 1. Initialize the ensemble with random states in [0, S-1]
    // ---------------------------------------
    std::vector<int> ensemble(number_of_sims, 0);
    pcg32 global_rng = get_random_number_generator();
    for (int j = 0; j < number_of_sims; j++){
        // pick a random integer in [0, S-1]
        ensemble[j] = global_rng() % S;
        //ensemble[j] = 0;
    }

    // ---------------------------------------
    // 2. Discrete time evolution
    // ---------------------------------------
    // choose a time step dt small enough so that (Tplus+Tminus)*dt << 1
    double dt = 0.01;  
    int num_steps = static_cast<int>(T / dt);

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    std::thread progress_thread(display_progress, std::ref(progress),num_steps);

    // For each time step, update all trajectories in parallel.
    // (For proper thread safety, one should give each thread its own RNG instance.)
    for (int step = 0; step < num_steps; step++){

        // Update trajectories in parallel
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < number_of_sims; j++){

            // If by chance this trajectory is absorbed, skip update (it will be cloned)
            if (ensemble[j] == S) continue;
            double state = ensemble[j];

            // Calculate rates (using the current state)
            double Tplus = calculate_T_plus_nonlinear_vm(N, state, Z, q);
            double Tminus = calculate_T_minus_nonlinear_vm(N, state, Z, q);
            double total_rate = Tplus + Tminus;

            // Probability that an event occurs in time dt
            double p_event = total_rate * dt;
            double r = draw_random_uniform(global_rng);
            if (r < p_event){

                // An event occurs: decide which event happens
                double r2 = draw_random_uniform(global_rng);
                if (r2 < Tplus / total_rate) state += 1;
                else state -= 1;

                // Ensure state remains in bounds:
                if (state < 0) state = 0;
                if (state > S) state = S;  // state==S means absorption.
                ensemble[j] = static_cast<int>(state);

            }
        } 

        // ---------------------------------------
        // 3. Cloning step: Replace any absorbed simulation by cloning another.
        // ---------------------------------------
        // Do this serially.
        for (int j = 0; j < number_of_sims; j++){
            if (ensemble[j] == S) {

                // Try to find a candidate that is not absorbed.
                int clone_index = -1;
                for (int attempt = 0; attempt < 10; attempt++){
                    int candidate = global_rng() % number_of_sims;
                    if (candidate == j) continue;
                    if (ensemble[candidate] != S){
                        clone_index = candidate;
                        break;
                    }
                }
                if (clone_index != -1){
                    ensemble[j] = ensemble[clone_index];
                }
                // If no candidate was found, leave ensemble[j] as S.
            }
        }

        progress++;

    }

    // Join the progress thread
    progress_thread.join();

    // ---------------------------------------
    // 4. Record the final states (only nonabsorbed ones) as the QSD sample.
    // ---------------------------------------
    std::unordered_map<int, int> state_counts;
    for (const auto& state : ensemble){
        if (state != S) {  // only count nonabsorbed trajectories
            state_counts[state]++;
        }
    }
    
    // count survived sims
    int survived_sims = 0;
    for (const auto& kv : state_counts) survived_sims += kv.second;
    
    std::vector<double> quasi_stationary_dist(S, 0.0);
    for (int i = 0; i < S; i++){
        if (state_counts.find(i) != state_counts.end()){
            quasi_stationary_dist[i] = static_cast<double>(state_counts[i]) / survived_sims;
        }
    }
    
    // Write the QSD to file
    write_single_vector_to_file(quasi_stationary_dist, folder_names, "q_sim_nonlinear_vm");
}

/**
 * @brief Solve analytically for the quasi-stationary state using recursive
 * techniques
 * 
 * @param q quasi-stationary distribution initialised to a uniform dist
 * @param all_Tplus All T+ values
 * @param all_Tminus All T- values
 * @param max_iter Maximum number of iterations
 * @param tolerance Tolerance for convergeace
 * @param size Number of susceptibles N-Z
 * @param verbose (default = true) print if needed
 */
void q_recursive_nonlinear_vm(std::vector<double> (&q), 
    const std::vector<double> (&all_Tplus), 
    const std::vector<double> (&all_Tminus), 
    const int& max_iter, const double& tolerance, const int& size,
    const bool& verbose = true)
{
    int iter{0};

    // this needs to be tuned
    double omega{0.5};

    // pre-calculate all alpha and beta coefficients
    std::vector<double> alpha(size), beta(size);
    for (int n = 1; n <= size-1; n++){
        alpha[n] = all_Tplus[n-1]/all_Tminus[n];
        beta[n] = all_Tplus[size-1]/all_Tminus[n];
    }

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    std::thread progress_thread(display_progress, std::ref(progress), max_iter);

    // repeatedly iterate until convergance
    std::vector<double> new_q(size, 0);
    std::vector<double> q_cum_sum(size, 0);
    double diff;
    for(int iter{0}; iter < max_iter; iter++)
    {

        // cummulative sum q
        std::partial_sum(q.begin(), q.end(), q_cum_sum.begin());

        // iterate
        new_q[0] = (1.-omega)*q[0] + omega*(all_Tminus[1]*q[1])/(all_Tplus[0] - all_Tplus[size-1]*q[size-1]);
        for (int n = 1; n <= size-1; n++){
            new_q[n] = (1.-omega)*q[n] + omega*(alpha[n]*q[n-1] - beta[n]*q[size-1]*q_cum_sum[n-1]);
        }

        // re-normalise
        double norm = 0;
        for (const auto& ele : new_q) norm += ele;
        for (auto& ele : new_q) ele /= norm;

        // calculate difference
        diff = 0;
        for (int i = 0; i <= size - 1; i++) diff += std::abs(new_q[i] - q[i]);

        // update q
        q = new_q;

        // break if converged
        if (diff < tolerance) break;

        progress++;
    }
    progress_thread.join();

    // print if requested
    if (verbose){
        if (iter < max_iter) std::cout << "Recurssive Converged!" << std::endl;
        else std::cout << "Recursive failed to converge :(" << std::endl;
    }

}

/**
 * @brief Solve for quasi-stationary distribution for the nonlinear vm
 *  using the Nifty approximation
 * 
 * @param N Population size
 * @return std::vector<double> quasidist
 */
std::vector<double> q_nifty_nonlinear_vm(
    const std::vector<double> (&all_Tplus), 
    const std::vector<double> (&all_Tminus), const int& size)
{
    std::vector<double> q_nifty(size, 0.);

    for (int n{0}; n <= size-1; n++){
        double temp{1};
        for (int j{0}; j <= n-1; j++)
            temp *= all_Tplus[j]/all_Tminus[j+1];
        q_nifty[n] = temp;
    }

    // calculate Q0
    double Q0{0};
    for (int n{0}; n <= size-1; n++){
        double temp{1};
        for (int j{0}; j <= n-1; j++)
            temp *= all_Tplus[j]/all_Tminus[j+1];
        Q0 += temp;
    }

    // divide all elements by Q0
    for (int n{0}; n <= size-1; n++) q_nifty[n] /= Q0;

    return q_nifty;
}


/**
 * @brief Calculate quasi-stationary distribution for the nonlinear voter model
 * with zealots in finite populations using recursive solve and the Nifty 
 * approx
 * 
 * @param q Param 
 * @param N Population size
 * @param Z Number of zealots
 */
void quasi_stationary_dist_analytic_nonlinear_vm(
    const double& q, const int& N, const int& Z)
{
    const int size = N-Z;

    // params
    std::stringstream qss;
    qss << "q=" << q;
    std::stringstream Zss;
    Zss << "Z=" << Z;

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "quasi_stationary_distribution", qss.str(), Zss.str()};

    // initialise quasi distribution as flat
    std::vector<double> q_recursive(size, 0.);
    for (int i{0}; i < size; i++) q_recursive[i] = 1./(size);

    // pre calculate all variables
    std::vector<double> all_Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);
    std::vector<double> all_Tminus = calculate_all_Tminus_nonlinear_vm(N, Z, q);

    // Perform iterative sole
    int max_iter = 1e4;
    double tolerance = 1e-10;
    q_recursive_nonlinear_vm(q_recursive, all_Tplus, all_Tminus, max_iter, 
        tolerance, size);

    // Perform nifty approx
    std::vector<double> q_nifty = q_nifty_nonlinear_vm(
        all_Tplus, all_Tminus, N-Z);

    // write data to file
    write_single_vector_to_file(q_recursive, folder_names, 
        "q_recursive_nonlinear_vm");
    write_single_vector_to_file(q_nifty, folder_names, 
        "q_nifty_nonlinear_vm");

}

/**
 * @brief Calculate the quasi stationary distribution numerically
 * for the nonlinear voter model in finite populations with absorbing zealots
 * 
 */
void quasi_stationary_dist_numeric_nonlinear_vm(
    const double& q, const int& N, const int& Z)
{
    // params
    std::stringstream qss;
    qss << "q=" << q;
    std::stringstream Zss;
    Zss << "Z=" << Z;

    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "quasi_stationary_distribution", qss.str(), Zss.str()};

    // initialise quasi distribution
    const int size = N-Z;
    Eigen::VectorXd quasi_dist = Eigen::VectorXd::Ones(size);
    quasi_dist /= quasi_dist.sum(); // Normalize to ensure it's a probability distribution

    // pre calculate all variables
    std::vector<double> all_Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);
    std::vector<double> all_Tminus = calculate_all_Tminus_nonlinear_vm(N, Z, q);
    Eigen::MatrixXd A = calculate_A_matrix(all_Tplus, all_Tminus, N, Z);

    double delta = 1e-5; // Time step
    int max_iter = 1e8;
    double tolerance = 1e-12;

    // Perform numerical integration
    explicit_euler(A, quasi_dist, all_Tplus[size-1], delta, max_iter, tolerance);

    // convert back to vector for writing
    std::vector<double> res;
    for (const auto& ele : quasi_dist) res.push_back(ele);

    // write data to file
    write_single_vector_to_file(res, folder_names, "q_masterly_nonlinear_vm");

}

/**
 * @brief For a variety of q values determine the Z value at which the 
 * recursive solution to the quasi stationary distribution is peaked at the
 * upper boundary
 * 
 * @param N system size
 */
void quasi_critical_z_nonlinear_vm(const int& N)
{
    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "quasi_stationary_distribution"};

    // containers
    std::vector<double> q_values = linspace(0.01, 3.01, 30); 
    std::vector<double> q_store; 
    std::vector<double> corresponding_Z; 

    // numerical integration params
    double delta = 0.001;
    int max_iter = 1e5;
    double tolerance = 1e-8;

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = q_values.size();

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // loop over q values
    for (const auto& q : q_values) {

        bool found = false;

        // loop over Z in [1, N-1] (i.e. always both susceptible and zealot)
        for (int Z = 1; Z <= N-1; Z++) {

            // initailsie flat distribution
            int size = N - Z;
            Eigen::VectorXd q_recursive = Eigen::VectorXd::Ones(size);
            q_recursive /= q_recursive.sum();

            // pre calculate all variables
            std::vector<double> all_Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);
            std::vector<double> all_Tminus = calculate_all_Tminus_nonlinear_vm(N, Z, q);
            Eigen::MatrixXd A = calculate_A_matrix(all_Tplus, all_Tminus, N, Z);
        
            // Perform numerical integration
            explicit_euler(A, q_recursive, all_Tplus[size-1], delta, max_iter, tolerance);

            // Determine the index of the maximum element.
            auto max_it = std::max_element(q_recursive.begin(), q_recursive.end());
            int max_index = std::distance(q_recursive.begin(), max_it);

            // Check if the maximum is at the upper boundary.
            if (max_index == size - 1) {
                q_store.push_back(q);
                corresponding_Z.push_back(Z);
                found = true;
                break;  // exit the Z loop once the condition is met.
            }
        }

        // Record -1 if no Z produces upper peal
        if (!found) {
            q_store.push_back(q);
            corresponding_Z.push_back(-1);  // or some flag value indicating no solution.
        }

        progress++;

    }
    progress_thread.join();

    // Write to file
    std::vector<std::vector<double>> res;
    for (size_t i{0}; i < q_store.size(); i++)
        res.push_back({q_store[i], corresponding_Z[i]});

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "quasi_critical_Z_varying_q");

}

/**
 * @brief Helper to return the quasi-stationary distribution, calculated
 * numerically, for the nonlinear vm with absorbing zealots
 * 
 * @param q nonlinear vm param
 * @param N system size
 * @param Z number of zealots
 * @return std::vector<double> Q_n 
 */
std::vector<double> quasi_helper(const double& q, const int& N, const int& Z)
{

    // initialise quasi distribution
    const int size = N-Z;
    Eigen::VectorXd quasi_dist = Eigen::VectorXd::Ones(size);
    quasi_dist /= quasi_dist.sum(); // Normalize to ensure it's a probability distribution

    // pre calculate all variables
    std::vector<double> all_Tplus = calculate_all_Tplus_nonlinear_vm(N, Z, q);
    std::vector<double> all_Tminus = calculate_all_Tminus_nonlinear_vm(N, Z, q);
    Eigen::MatrixXd A = calculate_A_matrix(all_Tplus, all_Tminus, N, Z);

    double delta, tolerance;
    int max_iter;
    if (q > 1 && Z >= 30){
        delta = 1e-5;
        max_iter = 1e6;
        tolerance = 1e-8;
    } else{
        delta = 1e-5; // Time step
        max_iter = 1e8;
        tolerance = 1e-12;
    }

    // Perform numerical integration
    explicit_euler(A, quasi_dist, all_Tplus[size-1], delta, max_iter, tolerance);

    // convert back to vector for writing
    std::vector<double> res;
    for (const auto& ele : quasi_dist) res.push_back(ele);

    return res;
}

/**
 * @brief For the nonlinear vm with absorbing zealots and in a finite 
 * population of size N determine the positions of the peak of the 
 * quasi-stationary distribution as we vary Z
 * 
 * @param q nonlinear vm param
 * @param N system size
 */
void quasi_stat_dist_peak_numeric_nonlinear_vm(const double& q, 
    const int& N)
{
    // set up folders for writing data
    std::vector<std::string> folder_names = {
        "nonlinear_vm", "quasi_stationary_distribution"};
    std::stringstream ss;
    ss << "q=" << q;
    folder_names.push_back(ss.str());

    // create z list
    const std::vector<int> z_values = arange(1, N-1);

    // container for results
    std::vector<std::vector<double>> res(N-1, std::vector<double>(2));

    // Initialize atomic progress counter
    std::atomic<int> progress(0);
    int total_iterations = z_values.size() - 1;

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), total_iterations);

    // Parallelize the loop using OpenMP
    // #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < z_values.size() - 1; i++) {
        int Z = z_values[i];

        // calculate the quasi-stationary distribution numerically
        std::vector<double> Q = quasi_helper(q, N, Z);

        // find the peak of the distribution - use std::max_element
        auto max_iter = std::max_element(Q.begin(), Q.end());
        double max_idx = std::distance(Q.begin(), max_iter);

        // store results
        res[i] = {static_cast<double>(Z)/N, max_idx/(N-Z-1.)};

        // Update progress counter - fixed atomic operation
        #pragma omp critical
        {
            progress++;
        }
    }
    progress_thread.join();

    // set final to avoid 0/0
    res[N-2] = {(N-1.)/N, 1.};

    // write data to file
    write_vectors_to_file_rowwise(res, folder_names, "quasi_peak_varying_Z");

}