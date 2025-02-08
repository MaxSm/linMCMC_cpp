# include "lin_ode_mcmc.h"

// Main program
int main(){

    // define parameters
    int iterations = 100000; //500000;
    double delta = 0.5;
    int steps = 2;
    double mu = 0.0202;
    double lambda = 0.02;

    // build pivotal matrices
    std::vector<double> linspace = make_linspace(delta, steps);
    key_matrices matrices = make_matrices(mu, lambda);

    // calculate
    results Results = calculate_results(iterations, matrices, linspace);
    
    // print out the result
    print_results(Results);

    // std::cout << "Hello world! \n" ; //debug point
    return 0;
}

