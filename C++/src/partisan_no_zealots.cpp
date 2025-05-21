#include "../includes/base_includes.hpp"
#include "../includes/utility.hpp"
#include "../includes/random.hpp"
#include "../includes/progress_bar.hpp"
#include "../includes/file_reading_and_writing.hpp"
#include "../includes/quasi_stationary_dist.hpp"

/**
 * @brief Calculate Wplus for all possible j values where j=N(delta+1/2) 
 * in the partisan vm without zealots
 * 
 * @param N population size
 * @param epsilon 
 * @return std::vector<double> all W+ 
 */
std::vector<double> calculate_all_Wplus_partisan_no_zealots(
    const double& N, const double& epsilon)
{
    std::vector<double> Wplus(N+1);
    for (int j{0}; j <= N; j++)
        Wplus[j] = (1.-static_cast<double>(j)/N)*j*(1.-2.*std::pow(epsilon, 2)*static_cast<double>(j)/N);

    return Wplus;
}

/**
 * @brief Calculate Wminus for all possible j values where j=N(delta+1/2) 
 * in the partisan vm without zealots
 * 
 * @param N population size
 * @param epsilon 
 * @return std::vector<double> all W-
 */
std::vector<double> calculate_all_Wminus_partisan_no_zealots(
    const double& N, const double& epsilon)
{
    std::vector<double> Wminus(N+1);
    for (int j{0}; j <= N; j++)
        Wminus[j] = (1.-static_cast<double>(j)/N)*j*(1.-2.*std::pow(epsilon, 2)*(1.-static_cast<double>(j)/N));

    return Wminus;
}

/**
 * @brief Solve for quasi-stationary distribution of partisan voter model
 * with no zealots (small epsilon) using the exact recursive method
 * 
 * @param max_iter Maximum number of iterations to perform
 * @param tolerance Tolerance level to stop iterating
 * @param N Size of systen
 */
void recursive_partisan_no_zealots(std::vector<double> (&q), 
    const std::vector<double> (&all_Wplus), 
    const std::vector<double> (&all_Wminus), 
    const int& max_iter, const double& tolerance, const int& N)
{
    int iter{0};

    // pre-calculate all alpha and beta coefficients
    std::vector<double> alpha(N), beta(N), gamma(N);
    for (int n = 1; n <= N-1; n++){
        alpha[n] = all_Wplus[n-1]/all_Wminus[n];
        beta[n] = all_Wminus[1]/all_Wminus[n];
        gamma[n] = all_Wplus[N-1]/all_Wminus[n];
    }

    // repeatedly iterate until convergance
    std::vector<double> new_q(N+1, 0);
    std::vector<double> q_cum_sum(N+1, 0);
    double diff;
    do {

        // cummulative sum q
        std::partial_sum(q.begin(), q.end(), q_cum_sum.begin());

        // iterate
        for (int n = 1; n <= N-1; n++){
            new_q[n] = alpha[n]*q[n-1]-gamma[n]*q[N-1]+(beta[n]*q[1]+gamma[n]*q[N-1])*(1-q_cum_sum[n-1]);
        }

        // re-normalise
        double norm = 0;
        for (int i{1}; i <= N-1; i++) norm += new_q[i];
        for (int i{1}; i <= N-1; i++) new_q[i] /= norm;

        // calculate difference
        diff = 0;
        for (int i = 0; i < N; i++) diff += std::abs(new_q[i] - q[i]);

        // update q
        q = new_q;

        iter++;

    } while (diff > tolerance && iter < max_iter);

    if (iter < max_iter) std::cout << "Recurssive Converged!" << std::endl;
    else std::cout << "Recursive failed to converge :(" << std::endl;

}

/**
 * @brief Solve for quasi-stationary distribution for the partisan voter
 * model with no zealots (small epsilon) using the Nifty approximation
 * 
 * @param N Population size
 * @return std::vector<double> quasidist
 */
std::vector<double> nifty_partisan_no_zealots(
    const std::vector<double> (&all_Wplus), 
    const std::vector<double> (&all_Wminus), const int& N)
{
    std::vector<double> q_nifty(N+1, 0.);

    for (int n{1}; n <= N-1; n++){
        double temp{1};
        for (int j{1}; j <= n-1; j++)
            temp *= all_Wplus[j]/all_Wminus[j+1];
        q_nifty[n] = temp;
    }

    // re-normalise
    double norm{0};
    for (int n{1}; n <= N-1; n++) norm += q_nifty[n];
    for (int n{1}; n <= N-1; n++) q_nifty[n] /= norm;

    return q_nifty;
}

/**
 * @brief Solve for the quasi-stationary distribution of the partisan voter
 * model without zealots (small epsilon) using both recursive solve and the 
 * nifty approx. Results are written to file.
 * 
 * @param epsilon Noise param (must be small) 
 * @param N populaiton size
 */
void quasi_stationary_dist_analytic_partisan_no_zealots(
    const double& epsilon, const double& N)
{
    // set up folders for writing data
    std::stringstream epsilon_ss;
    epsilon_ss << "epsilon=" << epsilon;
    std::vector<std::string> folder_names = {
        "partisan_no_zealots", "quasi_stationary_distribution", epsilon_ss.str()};

    // initialise quasi distribution as flat
    std::vector<double> q_recursive(N+1., 0.);
    for (int i{1}; i <= N-1; i++) q_recursive[i] = 1./(N-1.);

    // pre calculate all variables
    std::vector<double> all_Wplus = calculate_all_Wplus_partisan_no_zealots(
        N, epsilon);
    std::vector<double> all_Wminus = calculate_all_Wminus_partisan_no_zealots(
        N, epsilon);

    // Perform recursive solve
    int max_iter = 1e4;
    double tolerance = 1e-10;
    recursive_partisan_no_zealots(
        q_recursive, all_Wplus, all_Wminus, max_iter, tolerance, N);

    // Perform nifty approx
    std::vector<double> q_nifty = nifty_partisan_no_zealots(
        all_Wplus, all_Wminus, N);

    // write data to file
    write_single_vector_to_file(q_recursive, folder_names, 
        "q_recursive_partisan_no_zealots");
    write_single_vector_to_file(q_nifty, folder_names, 
        "q_nifty_partisan_no_zealots");


}















// /**
//  * @brief Form the transition matrix A an an Eigen matrix object for use in 
//  * numerically solving for the quasi stationary distribution
//  * 
//  * @param all_Wplus vector of all W+ values
//  * @param all_Wminus vector of all W- values
//  * @param N population size
//  * @return Eigen::MatrixXd 
//  */
// Eigen::MatrixXd calculate_A_matrix_partisan_no_zealots(
//     const std::vector<double>& all_Wplus, 
//     const std::vector<double>& all_Wminus, 
//     const int& N)
// {

//     // Initialize an Eigen matrix with zeros
//     Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N+1, N+1);

//     // Set diagonal elements
//     for (int i = 0; i <= N; i++) {
//         A(i, i) = -(all_Wplus[i] + all_Wminus[i]);
//     }

//     // Set upper diagonal elements
//     for (int i = 0; i <= N-1; i++) {
//         A(i, i + 1) = all_Wminus[i+1];
//     }

//     // Set lower diagonal elements
//     for (int i = 1; i <= N; i++) {
//         A(i, i - 1) = all_Wplus[i-1];
//     }

//     return A;
// }


// /**
//  * @brief Calculate the quasi stationary distribution numerically
//  * for the partisan voter model without zealots
//  * 
//  */
// void quasi_stationary_dist_numeric_partisan_no_zealots(
//     const double& epsilon, const int& N)
// {
//     // set up folders for writing data
//     std::stringstream epsilon_ss;
//     epsilon_ss << "epsilon=" << epsilon;
//     std::vector<std::string> folder_names = {
//         "partisan_no_zealots", "quasi_stationary_distribution", 
//         epsilon_ss.str()};

//     // initialise quasi distribution
//     Eigen::VectorXd quasi_dist = Eigen::VectorXd::Ones(N+1);
//     quasi_dist /= quasi_dist.sum(); // Normalize 

//     // pre calculate all variables
//     std::vector<double> all_Wplus = calculate_all_Wplus_partisan_no_zealots(
//         N, epsilon);
//     std::vector<double> all_Wminus = calculate_all_Wminus_partisan_no_zealots(
//         N, epsilon);
//     Eigen::MatrixXd A = calculate_A_matrix_partisan_no_zealots(
//        all_Wplus, all_Wminus, N);

//     double delta = 0.001; // Time step
//     int max_iter = 1e5;
//     double tolerance = 1e-8;

//     // Perform numerical integration
//     explicit_euler(A, quasi_dist, delta, max_iter, tolerance);

//     // convert back to vector for writing
//     std::vector<double> res;
//     for (const auto& ele : quasi_dist) res.push_back(ele);

//     // write data to file
//     write_single_vector_to_file(res, folder_names, "q_masterly_partisan_no_zealots");

// }





/**
 * @brief Determine quasi-stationary distribution for partisan voter model
 * without zealots using simulation techniques
 * 
 * @param epsilon Partisan parameter
 * @param N population size
 */
void quasi_stationary_dist_simulation_partisan_no_zealots(
    const double& epsilon, const int& N)
{

    // sim params
    const int number_of_sims = 1e5;
    const double T = 1000;
    const int j = N/2;  // related to delta via j=N(delta+1/2)

    // set up folders for writing data
    std::stringstream ss;
    ss << "epsilon=" << epsilon;
    std::vector<std::string> folder_names = {
        "partisan_no_zealots", "quasi_stationary_distribution", ss.str()};

    // pre calculate all variables
    const std::vector<double> all_Wplus = calculate_all_Wplus_partisan_no_zealots(
        N, epsilon);
    const std::vector<double> all_Wminus = calculate_all_Wminus_partisan_no_zealots(
        N, epsilon);

    // store container
    std::vector<double> j_store;

    // random number generator
    pcg32 source = get_random_number_generator();

    // Initialize atomic progress counter
    std::atomic<int> progress(0);

    // Start progress display thread
    std::thread progress_thread(display_progress, std::ref(progress), number_of_sims);

    // loop over all sims
    #pragma omp parallel for
    for (int i = 0; i < number_of_sims; i++){

        // initialise system
        int j_ind = j;
        double Wplus = all_Wplus[j_ind];
        double Wminus = all_Wminus[j_ind];
        double total_rate = Wplus + Wminus;

        // run for time T
        double t{0};
        double r1, delta;
        while (t < T){

            // determine update time
            r1 = draw_random_uniform(source);
            delta = -1.*log(r1)/total_rate;
            t += delta;

            // increase or decrease
            if (draw_random_uniform(source) <= Wplus/total_rate) j_ind += 1;
            else j_ind -= 1;

            // recalculate rates
            Wplus = all_Wplus[j_ind];
            Wminus = all_Wminus[j_ind];
            total_rate = Wplus + Wminus;

            // absorbing state
            if (j_ind == 0 || j_ind == N) break;

        }

        // store results
        if (j_ind != 0 && j_ind != N){
            #pragma omp critical
            {
                j_store.push_back(static_cast<int>(j_ind));
            }
        }

        progress++;
    }
    progress_thread.join();

    // count frequency each x appears
    std::unordered_map<int, int> mp;
    for (const auto& x : j_store) mp[x]++;

    // determine probaility distribution
    int survived_sims = j_store.size();
    std::vector<double> quasi_stationary_dist(N+1, 0.);
    for (int i{1}; i <= N-1; i++){

        // if number appeared
        if (mp.find(i) != mp.end()) quasi_stationary_dist[i] = static_cast<double>(mp[i])/survived_sims;
    }

    // write data to file
    write_single_vector_to_file(quasi_stationary_dist, folder_names, 
        "q_sim_partisan_no_zealots");

}