#include "../includes/base_includes.hpp"
#include "../includes/utility.hpp"
#include "../includes/random.hpp"
#include "../includes/file_reading_and_writing.hpp"
#include "nonlinear_vm.cpp"
#include "games.cpp"
#include "partisan.cpp"
#include "partisan_no_zealots.cpp"


int main()
{

    int Z = 24;
    int N = 100;

    // Nonlinear voter model
    double q = 3;
    //evolve_x_in_time_single_realisation_nonlinear_vm();
    //fixation_times_against_z_single_realisation_nonlinear_vm();
    //average_fixation_times_against_z_nonlinear_vm();
    //fixation_time_analytic_against_z_nonlinear_vm();
    //quasi_stationary_distribution_simulation_nonlinear_vm(q, N, Z);
    //quasi_stationary_dist_analytic_nonlinear_vm(q, N, Z);
    //quasi_stationary_dist_numeric_nonlinear_vm(q, N, Z);
    //quasi_critical_z_nonlinear_vm(N);
    quasi_stat_dist_peak_numeric_nonlinear_vm(q, N);

    // Evolutionary games
    // double alpha = 0.4;
    // average_fixation_times_against_z_games(N, alpha);
    // fixation_time_analytic_against_z_games(N, alpha);
    //quasi_critical_z_games(N);

    // Partisan voter model
    // double epsilon = 0.5;
    // average_fixation_times_against_z_partisan(epsilon, N);

    // Partisan vm no zealots (an in Llabres)
    // quasi_stationary_dist_analytic_partisan_no_zealots(epsilon, N);
    // quasi_stationary_dist_numeric_partisan_no_zealots(epsilon, N);
    // quasi_stationary_dist_simulation_partisan_no_zealots(epsilon, N);



    return 0;
}