# include "lin_ode_mcmc.h"

std::vector<double> make_linspace(double delta, int steps){
    std::vector<double> linspace;
    double tick = delta;
    for (int i = 0; i < steps; i++){
    linspace.push_back(tick);
    tick += delta;
    }
    return linspace;
};

const std::vector<std::vector<double>> A_matrix(double mu, double lambda){
    
    const std::vector<std::vector<double>> A_matrix = {
    {-10 * lambda, mu / 7, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {10 * lambda, -(9 * lambda + mu / 7), 2 * mu / 7, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 9 * lambda, -(8 * lambda + 2 * mu / 7), 3 * mu / 7, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 8 * lambda, -(7 * lambda + 3 * mu / 7), 4 * mu / 7, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 7 * lambda, -(6 * lambda + 4 * mu / 7), 5 * mu / 7, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 6 * lambda, -(5 * lambda + 5 * mu / 7), 6 * mu / 7, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 5 * lambda, -(4 * lambda + 6 * mu / 7), mu, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 4 * lambda, -(3 * lambda + mu), mu, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 3 * lambda, -(2 * lambda + mu), mu, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 2 * lambda, -(lambda + mu), mu},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, lambda, -mu}
    };
       
    return A_matrix;
}

const std::vector<std::vector<double>> P_matrix(){

    const std::vector<std::vector<double>> P_matrix = {
    {0.1, 0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.1, 0.1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0.1, 0.45, 0.15, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0.15, 0.4, 0.15, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0.2, 0.4, 0.1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0.35, 0.25, 0.1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0.4, 0.2, 0.1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0.4, 0.2, 0.1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0.45, 0.15, 0.1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0.45, 0.15, 0.1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0.1}
    };

    return P_matrix;
};

const std::vector<double> make_G(const std::vector<std::vector<double>> P){
    std::vector<double> G;
    for (int i = 0; i < P.size(); i++) {
    double sum = 0;
    for (int j = 0; j < 11; j++) {
        sum += P[i][j];
    }
    G.push_back(1 - sum);
    }
    const std::vector<double>& G_matrix = G; // std::move? :)
    return G_matrix;
};

const std::vector<double> make_p_0(){

    const std::vector<double> p_0 = {{0.05}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.05}};

    return p_0;    
};
const std::vector<double> make_H(){

    const std::vector<double> H = {{0.05}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.1}, {0.05}};

    return H;
};
std::vector<double> make_y_0(){

    std::vector<double> y_0 = {{0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};

    return y_0;
};

key_matrices make_matrices(double mu, double lambda){
    key_matrices matrices;
    matrices.A = A_matrix(mu, lambda);
    matrices.P = P_matrix();
    matrices.G = make_G(matrices.P);
    matrices.p_0 = make_p_0();  
    matrices.H = make_H();   
    matrices.y_0 = make_y_0();

    return matrices;
};

int initial_state(const std::vector<double> initial_distribution){
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double r = dis(gen);
    double p = initial_distribution[0];
    int state_0;
    for (int i = 0; i < initial_distribution.size(); i++){
        if (r < p) {
            state_0 = i;
            break;
        }
        else {
            p += initial_distribution[i + 1];
        }
    }

    return state_0;
};

int next_state(const std::vector<std::vector<double>> P, int previous_state){
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);     
    double r = dis(gen);
    std::vector<double> p_vec = P[previous_state];
    double p = p_vec[0];
    int state_i = -1;
    for (int i = 0; i < p_vec.size(); i++){
        if (r < p) {
            state_i = i;
            break;
        }
        else {
            p += p_vec[i + 1];
        }
    }

    return state_i;
};

std::vector<int> make_trajectory(const std::vector<double> initial_distribution, const std::vector<std::vector<double>> P){
    std::vector<int> trajectory;
    trajectory.push_back(initial_state(initial_distribution));
    while (trajectory.back() >= 0){
        trajectory.push_back(next_state(P, trajectory.back()));
    }
    trajectory.pop_back();
    return trajectory;
};

estimates calculate_J(std::vector<int> trajectory, const std::vector<double> p_0, const std::vector<double> H, \
const std::vector<double> G, const std::vector<double> y_0, const std::vector<std::vector<double>> P, \
const std::vector<std::vector<double>> A, double dt){
    std::random_device rd;  
    std::mt19937 gen(rd());
    estimates J_est;
    double j = H[trajectory[0]] / p_0[trajectory[0]];
    double j_sq = pow(j,2);
    for (int i = 1; i < trajectory.size(); i++){
        std::uniform_real_distribution<> dis(0.0, dt);     
        double delta = dis(gen);
        j = j * A[trajectory[i-1]][trajectory[i]] * delta / P[trajectory[i-1]][trajectory[i]]; //need t_matrix
        j_sq = j_sq * pow(A[trajectory[i-1]][trajectory[i]] * delta / P[trajectory[i-1]][trajectory[i]],2);
    }
    j = j * y_0[trajectory.back()] / G[trajectory.back()];
    j_sq = j_sq * pow(y_0[trajectory.back()] / G[trajectory.back()], 2);
    J_est.j = j;
    J_est.j_sq = j_sq;

    return J_est;
};

estimates_vector calculate_J_vector(int iterations, const std::vector<double> p_0, const std::vector<double> H, \
const std::vector<double> G, const std::vector<double> y_0, const std::vector<std::vector<double>> P, \
const std::vector<std::vector<double>> A, double dt){
    estimates_vector J_vector;
    for (int i = 0; i < iterations; i++){
        std::vector<int> trajectory = make_trajectory(p_0, P);
        J_vector.initial_state.push_back(trajectory[0]);
        estimates E = calculate_J(trajectory, p_0, H, G, y_0, P, A, dt);
        J_vector.J.push_back(E.j);
        J_vector.J_sq.push_back(E.j_sq);
    }
    return J_vector;
};

result calculate_result(estimates_vector est){
    result Res;
    std::vector<double> J_res (11, 0);
    std::vector<double> J_sq_res (11, 0);
    std::vector<double> D_res (11, 0);
    std::vector<double> St_err (11, 0);
    std::vector<double> Precision (11, 0);
    std::vector<int> count_array(11, 0);
    for (int i = 0; i < est.J.size(); i++){
        J_res[est.initial_state[i]] += est.J[i];
        J_sq_res[est.initial_state[i]] += est.J_sq[i];
        count_array[est.initial_state[i]] += 1;
    }
    for (int i = 0; i < J_res.size(); i++) {
        if (count_array[i] != 0){
            J_res[i] = J_res[i] / count_array[i];
            J_sq_res[i] = J_sq_res[i] / count_array[i];
        }
        else {
            J_res[i] = J_res[i] / 1;
            J_sq_res[i] = J_sq_res[i] / 1;
        }
        D_res[i] = J_sq_res[i] - pow(J_res[i], 2);
        St_err[i] = pow(D_res[i], 0.5);
        Precision[i] = 3 * St_err[i] / pow(est.J.size(), 0.5);
    }
    Res.solution = J_res;
    Res.standard_error = St_err;
    Res.precision = Precision;
    return Res;
};

results calculate_results(int iterations, key_matrices matrices, std::vector<double> linspace){
    results Results;
    double delta = linspace[0];
    for (int i = 0; i < linspace.size(); i++){
        estimates_vector J_vector = calculate_J_vector(iterations, matrices.p_0, matrices.H, matrices.G,\
        matrices.y_0, matrices.P, matrices.A, delta);    
        result Result = calculate_result(J_vector);
        matrices.y_0 = Result.solution;
        Results.solutions.push_back(Result);
        Results.time_stamp.push_back(linspace[i]);
    }    
    return Results;
};

void print_results(results Results){
     for (int j = 0; j < Results.solutions.size(); j++){
        std::cout << "Result" << "\t\t" << "Precision" << "\t" << "t = " << Results.time_stamp[j] <<"\n";
        for (int i = 0; i < Results.solutions[j].solution.size(); i++){
            std::cout << Results.solutions[j].solution[i] << "\t" << Results.solutions[j].precision[i] << "\n";
        }
        std::cout << "\n";
    }
};
