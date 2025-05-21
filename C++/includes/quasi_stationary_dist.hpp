#pragma once

#ifndef QUASI_STATIONARY_DIST
#define QUASI_STATIONARY_DIST

#include "base_includes.hpp"

/**
 * @brief Form the transition matrix A an an Eigen matrix object for use in 
 * numerically solving for the quasi stationary distribution
 * 
 * @param all_Tplus vector of all T+ values
 * @param all_Tminus vector of all T- values
 * @param N population size
 * @param Z number of +1 zealots
 * @return Eigen::MatrixXd 
 */
Eigen::MatrixXd calculate_A_matrix(
    const std::vector<double>& all_Tplus, 
    const std::vector<double>& all_Tminus, 
    const int& N, const int& Z)
{
    // Size of the matrix
    int size = N - Z;

    // Initialize an Eigen matrix with zeros
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);

    // Set diagonal elements
    for (int i = 0; i < size; i++) {
        A(i, i) = -(all_Tplus[i] + all_Tminus[i]);
    }

    // Set upper diagonal elements
    for (int i = 0; i < size - 1; i++) {
        A(i, i + 1) = all_Tminus[i+1];
    }

    // Set lower diagonal elements
    for (int i = 1; i < size; i++) {
        A(i, i - 1) = all_Tplus[i-1];
    }

    return A;
}

/**
 * @brief Perform explict (forward) Euler numerical integration using 
 * the Eigen package
 * 
 * @param A Transition matrix with size (N-Z, N-Z)
 * @param q quasi-stationary row vector of size (N-Z)
 * @param Tplus single element T+[N-Z-1]
 * @param delta time step
 * @param max_iter max number of iterations
 * @param tolerance threshold for stopping numerical integration
 */
void explicit_euler(const Eigen::MatrixXd& A, Eigen::VectorXd& q, 
    const double& Tplus, const double& delta, const int& max_iter, 
    const double& tolerance)
{
    int size = q.size();
    Eigen::VectorXd q_next(size);
    double relative_change{1};
    for (int i = 0; i < max_iter; i++) {
        q_next.noalias() = q + delta * (A*q + Tplus*q[size-1]*q);
        relative_change = (q_next-q).cwiseAbs().sum();
        q = std::move(q_next);
        if (relative_change < tolerance) {
            std::cout << "Converged!" << std::endl;
            break;
        }
    }

    // check for convergence
    if (relative_change > tolerance){
        std::cout << "Warning: Not converged!" << std::endl;
    }
}




#endif