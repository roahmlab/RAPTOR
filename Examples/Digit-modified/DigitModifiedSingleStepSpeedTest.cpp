#include "DigitModifiedSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace RAPTOR;
using namespace DigitModified;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit-modified/data/";

const int degree_range[] = {5, 6, 7, 8};
const int N_range[] = {16, 20, 24, 28};
const std::string mu_strategy_str[] = {"adaptive", "monotone"};
const TimeDiscretization time_discret[] = {Uniform, Chebyshev};
const std::string time_discretization_str[] = {"Uniform", "Chebyshev"};

int main(int argc, char* argv[]) {
    // Define robot model
    const std::string urdf_filename = "../Robots/digit-v3-modified/digit-v3-modified.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    if (model.nq != 20) {
        throw std::invalid_argument("Error: Incorrect number of joints in the robot model!");
    }
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 3, 3, 3,
             3, 3, -3, 3, 3, 3, 3;
    
    // ignore all motor dynamics
    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();

    // load settings
    const double T = 0.4;
    int N = 16;
    int degree = 5;
    
    GaitParameters gp;
    gp.swingfoot_midstep_z_des = 0.30;
    gp.swingfoot_begin_y_des = 0.40;
    gp.swingfoot_end_y_des = -0.40;

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-5);
    app->Options()->SetNumericValue("max_wall_time", 200.0);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    std::ofstream experiment_output(filepath + "speed_test.txt");

    for (int degree_choice = 0; degree_choice < 4; degree_choice++) {
        for (int mu_strategy_choice = 0; mu_strategy_choice < 2; mu_strategy_choice++) {
            for (int discretization_choice = 0; discretization_choice < 2; discretization_choice++) {
                degree = degree_range[degree_choice];
                N = N_range[degree_choice];
                app->Options()->SetStringValue("mu_strategy", mu_strategy_str[mu_strategy_choice].c_str());
                TimeDiscretization time_discretization = time_discret[discretization_choice];

                std::cerr << "EXPERIMENT: Degree: " << degree << ", mu_strategy: " << mu_strategy_str[mu_strategy_choice] << ", time discretization: " << time_discretization_str[discretization_choice] << std::endl;
                // experiment_output << "EXPERIMENT: Degree: " << degree << ", mu_strategy: " << mu_strategy_str[mu_strategy_choice] << ", time discretization: " << time_discretization_str[discretization_choice] << std::endl;

                Eigen::VectorXd z = Utils::initializeEigenMatrixFromFile(filepath + "robustness_test_solution_" + std::to_string(degree) + ".txt");
                
                SmartPtr<DigitModifiedSingleStepOptimizer> mynlp = new DigitModifiedSingleStepOptimizer();
                try {
                    mynlp->set_parameters(z,
                                          T,
                                          N,
                                          time_discretization,
                                          degree,
                                          model,
                                          jtype,
                                          gp);
                    mynlp->constr_viol_tol = 1e-5;
                }
                catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
                }

                // Initialize the IpoptApplication and process the options
                ApplicationReturnStatus status;
                status = app->Initialize();
                if( status != Solve_Succeeded ) {
                    throw std::runtime_error("Error during initialization of optimization!");
                }

                // Solve the optimization problem
                double solve_time = 0.0;
                try {
                    auto start = std::chrono::high_resolution_clock::now();

                    status = app->OptimizeTNLP(mynlp);

                    auto end = std::chrono::high_resolution_clock::now();
                    solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;

                    experiment_output << " & & & " << mynlp->obj_value_copy << " & " << mynlp->final_constr_violation << " & ";
                }
                catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    throw std::runtime_error("Error solving optimization problem! Check previous error message!");
                }

                // Evaluate the solution on a finer time discretization
                try {
                    SmartPtr<DigitModifiedSingleStepOptimizer> testnlp = new DigitModifiedSingleStepOptimizer();
                    testnlp->set_parameters(z,
                                            T,
                                            2000,
                                            TimeDiscretization::Uniform,
                                            degree,
                                            model,
                                            jtype,
                                            gp);
                    Index n, m, nnz_jac_g, nnz_h_lag;
                    TNLP::IndexStyleEnum index_style;
                    testnlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
                    Number ztry[testnlp->numVars], x_l[testnlp->numVars], x_u[testnlp->numVars];
                    Number g[testnlp->numCons], g_lb[testnlp->numCons], g_ub[testnlp->numCons];
                    for (int i = 0; i < testnlp->numVars; i++) {
                        ztry[i] = mynlp->solution[i];
                    }
                    testnlp->get_bounds_info(testnlp->numVars, x_l, x_u, testnlp->numCons, g_lb, g_ub);
                    testnlp->eval_g(testnlp->numVars, ztry, false, testnlp->numCons, g);
                    testnlp->summarize_constraints(testnlp->numCons, g, true);
                    experiment_output << testnlp->final_constr_violation << " & " << solve_time << " \\\\" << std::endl;
                }
                catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    throw std::runtime_error("Error evaluating the solution on a finer time discretization! Check previous error message!");
                }
            }
        }
    }

    experiment_output.close();

    return 0;
}
