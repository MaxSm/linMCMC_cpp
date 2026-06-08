# ifndef LIN_ODE_MCMC_H
# define LIN_ODE_MCMC_H

# include <iostream>
# include <fstream>
# include <vector>
# include <array>
# include <random>
# include <cmath>


struct key_matrices{
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> P;
    std::vector<double> G;
    std::vector<double> p_0;  
    std::vector<double> H;   
    std::vector<double> y_0;
};

struct estimates{
    double j;
    double j_sq;
};

struct estimates_vector{
    std::vector<double> J;
    std::vector<double> J_sq;
    std::vector<int> initial_state;
};

struct result{
    std::vector<double> solution;
    std::vector<double> standard_error;
    std::vector<double> precision;
};

struct results{
    std::vector<result> solutions;
    std::vector<double> time_stamp;
};

std::vector<double> make_linspace(double delta, int steps);

const std::vector<std::vector<double>> A_matrix(double mu, double lambda);

const std::vector<std::vector<double>> P_matrix();

const std::vector<double> make_G(const std::vector<std::vector<double>> P);

const std::vector<double> make_p_0();

const std::vector<double> make_H();

std::vector<double> make_y_0();

key_matrices make_matrices(double mu, double lambda);

int initial_state(const std::vector<double> initial_distribution);

int next_state(const std::vector<std::vector<double>> P, int previous_state);

std::vector<int> make_trajectory(const std::vector<double> initial_distribution, const std::vector<std::vector<double>> P);

estimates calculate_J(std::vector<int> trajectory, const std::vector<double> p_0, const std::vector<double> H, \
const std::vector<double> G, const std::vector<double> y_0, const std::vector<std::vector<double>> P, \
const std::vector<std::vector<double>> A, double dt);

estimates_vector calculate_J_vector(int iterations, const std::vector<double> p_0, const std::vector<double> H, \
const std::vector<double> G, const std::vector<double> y_0, const std::vector<std::vector<double>> P, \
const std::vector<std::vector<double>> A, double dt);

result calculate_result(estimates_vector estimates);

results calculate_results(int iterations, key_matrices matrices, std::vector<double> linspace);

void print_results(results Results);

#endif