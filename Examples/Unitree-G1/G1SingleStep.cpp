#include "G1SingleStepOptimizer.h"

#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <yaml-cpp/yaml.h>
#include <iomanip>

using namespace RAPTOR;
using namespace G1;
using namespace Ipopt;

const std::string filepath = "../Examples/Unitree-G1/data/";

int main(int argc, char* argv[]) {
    // define robot model
    const std::string urdf_filename = "../Robots/unitree-g1/g1_22dof_floatingbase.urdf";
    // const std::string urdf_filename = "../Robots/unitree-g1/g1_12dof_floatingbase.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // lock the arm joints at a fixed position
    pinocchio::Model model_armfixed;
    std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id(10);
    for (int i = 0; i < 10; i++) {
        list_of_joints_to_lock_by_id[i] = (pinocchio::JointIndex)(i + 12 + 6 + 1);
    }
    Eigen::VectorXd q_armfixed = Eigen::VectorXd::Zero(model.nq);
    q_armfixed(12) = 0.6; // left_shoulder_pitch_joint
    q_armfixed(13) = 0.2; // left_shoulder_roll_joint
    q_armfixed(15) = 0.3; // left_elbow_pitch_joint
    q_armfixed(17) = 0.6; // right_shoulder_pitch_joint
    q_armfixed(18) = -0.2; // right_shoulder_roll_joint
    q_armfixed(20) = 0.3; // right_elbow_pitch_joint
    pinocchio::buildReducedModel(model, list_of_joints_to_lock_by_id, q_armfixed, model_armfixed);
    model = model_armfixed;
    
    // ignore all motor dynamics
    model.armature.setZero();
    model.damping.setZero();
    model.friction.setZero();

    // load settings
    YAML::Node config;

    const double T = 0.4;
    TimeDiscretization time_discretization = Uniform;
    int N = 14;
    int degree = 5;
    
    GaitParameters gp;

    try {
        config = YAML::LoadFile("../Examples/Unitree-G1/singlestep_optimization_settings.yaml");

        N = config["N"].as<int>();
        degree = config["degree"].as<int>();
        std::string time_discretization_str = config["time_discretization"].as<std::string>();
        time_discretization = (time_discretization_str == "Uniform") ? Uniform : Chebyshev;

        gp.swingfoot_midstep_z_des = config["swingfoot_midstep_z_des"].as<double>();
        gp.swingfoot_begin_x_des = config["swingfoot_begin_x_des"].as<double>();
        gp.swingfoot_end_x_des = config["swingfoot_end_x_des"].as<double>();
    } 
    catch (std::exception& e) {
        std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
    }
    
    if (argc > 1) {
        char* end = nullptr;
        std::srand((unsigned int)std::strtoul(argv[1], &end, 10));
    }
    else {
        std::srand(std::time(nullptr));
    }
    Eigen::VectorXd z = 0.1 * Eigen::VectorXd::Random((degree + 1) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS);
    
    SmartPtr<G1SingleStepOptimizer> mynlp = new G1SingleStepOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              time_discretization,
                              degree,
                              model,
                              gp);
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
        std::ofstream solution(filepath + "solution-g1-forward.txt");

        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }
        solution.close();

        std::ofstream trajectory(filepath + "trajectory-g1.txt");
        trajectory << std::setprecision(20);
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->cidPtr_->q(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->cidPtr_->v(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->cidPtr_->a(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->cidPtr_->tau(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->cidPtr_->lambda(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        trajectory.close();
    }

    return 0;
}
