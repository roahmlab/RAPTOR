#include "DigitMultipleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace IDTO;
using namespace Digit;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit/data/";

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

    const int NSteps = 2;
    const double T = 0.4;
    const TimeDiscretization time_discretization = Uniform;
    const int N = 8;
    const int degree = 5;
    // const std::string output_name = std::string(argv[1]) + "-" + std::string(argv[2]);
    const std::string output_name = "Two-8-5-Uniform";

    GaitParameters gp;
    gp.swingfoot_midstep_z_des = 0.15;
    gp.swingfoot_begin_y_des = 0.00;
    gp.swingfoot_end_y_des = -0.00;
    // gp.swingfoot_midstep_z_des = std::atof(argv[2]);
    // gp.swingfoot_begin_y_des = std::atof(argv[1]);
    // gp.swingfoot_end_y_des = - std::atof(argv[1]);
    
    Eigen::VectorXd z = Utils::initializeEigenMatrixFromFile(filepath + "initial-digit-Bezier-Two-14-5-Uniform.txt");

    SmartPtr<DigitMultipleStepOptimizer> mynlp = new DigitMultipleStepOptimizer();
    try {
	    mynlp->set_parameters(NSteps,
                              z,
                              T,
                              N,
                              time_discretization,
                              degree,
                              model,
                              jtype,
                              gp);
        mynlp->constr_viol_tol = 1e-4;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetNumericValue("max_wall_time", 1e-4);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("nlp_scaling_method", "none");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt_digit.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("point_perturbation_radius", 1e-3);
    // // app->Options()->SetIntegerValue("derivative_test_first_index", 168);
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-4);

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
    // if (mynlp->solution.size() == mynlp->numVars) {
        // std::ofstream solution(filepath + "solution-digit-Bezier-" + output_name + ".txt");
        // solution << std::setprecision(20);
        // for (int i = 0; i < mynlp->numVars; i++) {
        //     solution << mynlp->solution[i] << std::endl;
        // }
        // solution.close();

        std::ofstream trajectory(filepath + "trajectory-digit-Bezier-" + output_name + ".txt");
        trajectory << std::setprecision(20);
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < mynlp->cidPtr_->N; j++) {
                trajectory << mynlp->cidPtr_->q(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < mynlp->cidPtr_->N; j++) {
                trajectory << mynlp->cidPtr_->v(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < mynlp->cidPtr_->N; j++) {
                trajectory << mynlp->cidPtr_->a(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
            for (int j = 0; j < mynlp->cidPtr_->N; j++) {
                trajectory << mynlp->cidPtr_->tau(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
            for (int j = 0; j < mynlp->cidPtr_->N; j++) {
                trajectory << mynlp->cidPtr_->lambda(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        trajectory.close();
    // }

    return 0;
}
