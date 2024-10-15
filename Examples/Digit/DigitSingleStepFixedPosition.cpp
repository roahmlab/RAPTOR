#include "DigitSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <yaml-cpp/yaml.h>
#include <iomanip>

using namespace RAPTOR;
using namespace Digit;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit/data/";

int main(int argc, char* argv[]) {
    // define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;
    
    // ignore friction for now
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

    // the robot start from this initial condition:
    Eigen::VectorXd q0(NUM_INDEPENDENT_JOINTS);
    q0 << 0.35717257333012147890727305821201,
          0.0068511286518564409528386782710641,
          0.34271261983060330447159458344686,
          0.34477720686111551495756089025235,
          -0.13786708397138611181809153549693,
          0.12782183922354697913625898308965,
          -0.35247716987640814734206173852726,
          0.027851543514003956891755819924583,
          -0.34806008324315867996645579296455,
          -0.34836351024799677711385470502137,
          0.14676360138953425948926678756834,
          -0.12790655894451277707624115009821;
    Eigen::VectorXd q_d0(NUM_INDEPENDENT_JOINTS);
    q_d0.setZero();

    // load settings
    YAML::Node config;

    const double T = 0.35;
    TimeDiscretization time_discretization = Uniform;
    int N = 14;
    int degree = 5;
    
    GaitParameters gp;

    try {
        config = YAML::LoadFile("../Examples/Digit/singlestep_optimization_settings.yaml");

        N = config["N"].as<int>();
        degree = config["degree"].as<int>();
        std::string time_discretization_str = config["time_discretization"].as<std::string>();
        time_discretization = (time_discretization_str == "Uniform") ? Uniform : Chebyshev;

        gp.swingfoot_midstep_z_des = config["swingfoot_midstep_z_des"].as<double>();
        gp.swingfoot_begin_y_des = config["swingfoot_begin_y_des"].as<double>();
        gp.swingfoot_end_y_des = config["swingfoot_end_y_des"].as<double>();
    } 
    catch (std::exception& e) {
        std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
    }

    // const std::string output_name = std::string(argv[1]) + "-" + std::string(argv[2]);
    
    // Eigen::VectorXd z = Utils::initializeEigenMatrixFromFile(filepath + "initial-digit-Bezier-14-5-Uniform.txt");
    // if (argc > 1) {
    //     char* end = nullptr;
    //     std::srand((unsigned int)std::strtoul(argv[1], &end, 10));
    // }
    // else {
    //     std::srand(std::time(nullptr));
    // }
    // Eigen::VectorXd z = 0.2 * Eigen::VectorXd::Random((degree) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS).array() - 0.1;
    Eigen::VectorXd z = Eigen::VectorXd::Zero((degree) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS);
    
    SmartPtr<DigitSingleStepOptimizer> mynlp = new DigitSingleStepOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              time_discretization,
                              degree,
                              model,
                              gp,
                              'L',
                              Transform(3, -M_PI_2),
                              false,
                              q0);
        mynlp->constr_viol_tol = config["constr_viol_tol"].as<double>();
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    try {
        app->Options()->SetNumericValue("tol", config["tol"].as<double>());
        app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
        app->Options()->SetNumericValue("max_wall_time", config["max_wall_time"].as<double>());

        app->Options()->SetIntegerValue("max_iter", config["max_iter"].as<int>());
        // char* end = nullptr;
        // app->Options()->SetIntegerValue("max_iter", (unsigned int)std::strtoul(argv[1], &end, 10));

        app->Options()->SetNumericValue("obj_scaling_factor", config["obj_scaling_factor"].as<double>());
        app->Options()->SetIntegerValue("print_level", config["print_level"].as<double>());
        app->Options()->SetStringValue("mu_strategy", config["mu_strategy"].as<std::string>().c_str());
        app->Options()->SetStringValue("linear_solver", config["linear_solver"].as<std::string>().c_str());
        app->Options()->SetStringValue("ma57_automatic_scaling", "yes");

        if (mynlp->enable_hessian) {
            app->Options()->SetStringValue("hessian_approximation", "exact");
        }
        else {
            app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        }
        // app->Options()->SetStringValue("nlp_scaling_method", "none");
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error setting optimization options! Check previous error message!");
    }

    if (config["gredient_check"].as<bool>()) {
        app->Options()->SetStringValue("output_file", "ipopt.out");
        if (mynlp->enable_hessian) {
            app->Options()->SetStringValue("derivative_test", "second-order");
        }
        else {
            app->Options()->SetStringValue("derivative_test", "first-order");
        }
        app->Options()->SetNumericValue("point_perturbation_radius", 1e-3);
        // app->Options()->SetIntegerValue("derivative_test_first_index", 168);
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-4);
    }

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    // char numBuffer[10];
    // sprintf(numBuffer, "%.1f", abs(gp.swingfoot_end_y_des));
    // std::ofstream experiment_output(filepath + "speed_output_" + numBuffer + ".txt", std::ofstream::out | std::ofstream::app);

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        double solve_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6;

        std::cout << "Data needed for comparison: " << mynlp->obj_value_copy << ' ' << mynlp->final_constr_violation << ' ' << solve_time << std::endl;
        // experiment_output << mynlp->obj_value_copy << ' ' << mynlp->final_constr_violation << ' ' << solve_time << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        // Evaluate the solution on a finer time discretization
        try {
            const double dt_sim = 5e-4;
            const int N_simulate = T / dt_sim + 1;

            SmartPtr<DigitSingleStepOptimizer> testnlp = new DigitSingleStepOptimizer();
            testnlp->display_info = false;
            testnlp->set_parameters(z,
                                    T,
                                    N_simulate,
                                    TimeDiscretization::Uniform,
                                    degree,
                                    model,
                                    gp,
                                    'L',
                                    Transform(3, -M_PI_2),
                                    false,
                                    q0);
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
            testnlp->summarize_constraints(testnlp->numCons, g, false);

            char numBuffer[10];
            sprintf(numBuffer, "%.1f", abs(gp.swingfoot_end_y_des));
            std::ofstream solution(filepath + "solution-digit-forward-" + numBuffer + ".txt");
            for (int j = 0; j < N_simulate; j++) {
                for (int i = 0; i < NUM_JOINTS; i++) {
                    solution << testnlp->cidPtr_->q(j)(i) << ' ';
                }
                for (int i = 0; i < NUM_JOINTS; i++) {
                    solution << testnlp->cidPtr_->v(j)(i) << ' ';
                }
                for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
                    solution << testnlp->cidPtr_->tau(j)(i) << ' ';
                }
                solution << std::endl;
            }
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            throw std::runtime_error("Error evaluating the solution on a finer time discretization! Check previous error message!");
        }

        // std::ofstream solution(filepath + 
        //                        "robustness_test_solution_" + 
        //                        std::to_string(degree) + 
        //                        ".txt");

        // solution << std::setprecision(20);
        // for (int i = 0; i < mynlp->numVars; i++) {
        //     solution << mynlp->solution[i] << std::endl;
        // }
        // solution.close();

        // std::ofstream trajectory(filepath + "trajectory-digit-Bezier-" + output_name + ".txt");
        // trajectory << std::setprecision(20);
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->cidPtr_->q(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->cidPtr_->v(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->cidPtr_->a(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->cidPtr_->tau(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->cidPtr_->lambda(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // trajectory.close();
    }

    return 0;
}
