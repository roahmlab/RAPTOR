#include "DigitMultipleStepOptimizer.h"

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

    // manually import motor inertia 
    model.armature(model.getJointId("left_hip_roll") - 1) = 0.173823936;
    model.armature(model.getJointId("left_hip_yaw") - 1) = 0.067899975;
    model.armature(model.getJointId("left_hip_pitch") - 1) = 0.1204731904;
    model.armature(model.getJointId("left_knee") - 1) = 0.1204731904;
    model.armature(model.getJointId("left_toe_A") - 1) = 0.036089475;
    model.armature(model.getJointId("left_toe_B") - 1) = 0.036089475;
    model.armature(model.getJointId("right_hip_roll") - 1) = 0.173823936;
    model.armature(model.getJointId("right_hip_yaw") - 1) = 0.067899975;
    model.armature(model.getJointId("right_hip_pitch") - 1) = 0.1204731904;
    model.armature(model.getJointId("right_knee") - 1) = 0.1204731904;
    model.armature(model.getJointId("right_toe_A") - 1) = 0.036089475;
    model.armature(model.getJointId("right_toe_B") - 1) = 0.036089475;

    // load settings
    YAML::Node config;

    const double T = 0.4;
    int NSteps = 2;
    TimeDiscretization time_discretization = Uniform;
    int N = 14;
    int degree = 5;
    
    std::vector<GaitParameters> gps(NSteps);

    try {
        config = YAML::LoadFile("../Examples/Digit/multiplestep_optimization_settings.yaml");

        NSteps = config["NSteps"].as<int>();
        N = config["N"].as<int>();
        degree = config["degree"].as<int>();
        std::string time_discretization_str = config["time_discretization"].as<std::string>();
        time_discretization = (time_discretization_str == "Uniform") ? Uniform : Chebyshev;

        auto swingfoot_midstep_z_des = config["swingfoot_midstep_z_des"].as<std::vector<double>>();
        auto swingfoot_begin_y_des = config["swingfoot_begin_y_des"].as<std::vector<double>>();
        auto swingfoot_end_y_des = config["swingfoot_end_y_des"].as<std::vector<double>>();

        if (swingfoot_midstep_z_des.size() != NSteps || 
            swingfoot_begin_y_des.size() != NSteps || 
            swingfoot_end_y_des.size() != NSteps) {
            throw std::runtime_error("Error parsing YAML file: Incorrect number of gait parameters!");
        }

        gps.resize(NSteps);
        for (int i = 0; i < NSteps; i++) {
            gps[i].swingfoot_midstep_z_des = swingfoot_midstep_z_des[i];
            gps[i].swingfoot_begin_y_des = swingfoot_begin_y_des[i];
            gps[i].swingfoot_end_y_des = swingfoot_end_y_des[i];
        }
    } 
    catch (std::exception& e) {
        std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
        throw std::runtime_error("Error parsing YAML file! Check previous error message!");
    }
    
    Eigen::VectorXd z_onestep = Utils::initializeEigenMatrixFromFile(filepath + "initial-digit.txt");
    Eigen::VectorXd z0(z_onestep.size() * NSteps);
    for (int i = 0; i < NSteps; i++) {
        z0.segment(i * z_onestep.size(), z_onestep.size()) = z_onestep;
    }
    // add noise to initial guess to explore the solution space
    // std::srand(std::time(nullptr));
    // z = z + 0.01 * Eigen::VectorXd::Random(z.size());

    SmartPtr<DigitMultipleStepOptimizer> mynlp = new DigitMultipleStepOptimizer();
    try {
	    mynlp->set_parameters(NSteps,
                              z0,
                              T,
                              N,
                              time_discretization,
                              degree,
                              model,
                              gps,
                              true);
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

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        double solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        
        std::cout << "Data needed for comparison: " << mynlp->obj_value_copy << ' ' << mynlp->final_constr_violation << ' ' << solve_time << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        std::ofstream solution(filepath + "solution-digit-multiple-step.txt");
        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }
        solution.close();

        for (int p = 0; p < NSteps; p++) {
            std::ofstream trajectory(filepath + "trajectory-digit-multiple-step-" + std::to_string(p) + ".txt");
            trajectory << std::setprecision(20);
            const auto& cidPtr_ = mynlp->stepOptVec_[p]->cidPtr_;
            for (int i = 0; i < NUM_JOINTS; i++) {
                for (int j = 0; j < cidPtr_->N; j++) {
                    trajectory << cidPtr_->q(j)(i) << ' ';
                }
                trajectory << std::endl;
            }
            for (int i = 0; i < NUM_JOINTS; i++) {
                for (int j = 0; j < cidPtr_->N; j++) {
                    trajectory << cidPtr_->v(j)(i) << ' ';
                }
                trajectory << std::endl;
            }
            for (int i = 0; i < NUM_JOINTS; i++) {
                for (int j = 0; j < cidPtr_->N; j++) {
                    trajectory << cidPtr_->a(j)(i) << ' ';
                }
                trajectory << std::endl;
            }
            for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
                for (int j = 0; j < cidPtr_->N; j++) {
                    trajectory << cidPtr_->tau(j)(i) << ' ';
                }
                trajectory << std::endl;
            }
            for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
                for (int j = 0; j < cidPtr_->N; j++) {
                    trajectory << cidPtr_->lambda(j)(i) << ' ';
                }
                trajectory << std::endl;
            }
            trajectory.close();
        }
    }

    return 0;
}
