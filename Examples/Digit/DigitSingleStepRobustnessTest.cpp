#include "DigitSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace IDTO;
using namespace Digit;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit/data/robustness_test/";

const int degree_range[] = {5, 6, 7, 8};
const int N_range[] = {14, 18, 22, 26};
const std::string mu_strategy_str[] = {"adaptive", "monotone"};
const std::string linear_solver_str[] = {"ma27", "ma57", "ma86"};

int main(int argc, char* argv[]) {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);

    // define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.806;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3,
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3;
    
    // ignore friction for now
    model.friction.setZero();

    // manually import motor inertia 
    model.rotorInertia(model.getJointId("left_hip_roll") - 1) = 0.173823936;
    model.rotorInertia(model.getJointId("left_hip_yaw") - 1) = 0.067899975;
    model.rotorInertia(model.getJointId("left_hip_pitch") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("left_knee") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("left_toe_A") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("left_toe_B") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("right_hip_roll") - 1) = 0.173823936;
    model.rotorInertia(model.getJointId("right_hip_yaw") - 1) = 0.067899975;
    model.rotorInertia(model.getJointId("right_hip_pitch") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("right_knee") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("right_toe_A") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("right_toe_B") - 1) = 0.036089475;

    // load settings
    const double T = 0.4;
    const TimeDiscretization time_discretization = Chebyshev;
    int N = 16;
    int degree = 5;
    
    GaitParameters gp;

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-5);
    app->Options()->SetNumericValue("max_wall_time", 200.0);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    for (int degree_choice = 0; degree_choice < 4; degree_choice++) {
        degree = degree_range[degree_choice];
        N = N_range[degree_choice];
        for (int mu_strategy_choice = 0; mu_strategy_choice < 2; mu_strategy_choice++) {
            app->Options()->SetStringValue("mu_strategy", mu_strategy_str[mu_strategy_choice].c_str());
            for (int linear_solver_choice = 0; linear_solver_choice < 2; linear_solver_choice++) {
                app->Options()->SetStringValue("linear_solver", linear_solver_str[linear_solver_choice].c_str());

                std::cerr << "EXPERIMENT: Degree: " << degree << ", mu_strategy: " << mu_strategy_str[mu_strategy_choice] << ", linear_solver: " << linear_solver_str[linear_solver_choice] << std::endl;
                
                std::ofstream experiment_output(filepath + 
                                                "robustness_test_" + 
                                                std::to_string(degree) + 
                                                "_" +
                                                mu_strategy_str[mu_strategy_choice] +
                                                "_" +
                                                linear_solver_str[linear_solver_choice] +
                                                ".txt");

                for (int test_id = 1; test_id <= 100; test_id++) {
                    std::srand(test_id);
                    Eigen::VectorXd z = 0.2 * Eigen::VectorXd::Random((degree + 1) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS).array() - 0.1;
                    
                    SmartPtr<DigitSingleStepOptimizer> mynlp = new DigitSingleStepOptimizer();
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

                    try {
                        auto start = std::chrono::high_resolution_clock::now();

                        // Ask Ipopt to solve the problem
                        status = app->OptimizeTNLP(mynlp);

                        auto end = std::chrono::high_resolution_clock::now();
                        double solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;

                        experiment_output << mynlp->obj_value_copy << ' ' << mynlp->final_constr_violation << ' ' << solve_time << std::endl;
                    }
                    catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
                    }
                }

                experiment_output.close();
            }
        }
    }

    return 0;
}
