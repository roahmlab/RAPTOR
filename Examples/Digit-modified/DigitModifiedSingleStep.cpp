#include "DigitModifiedSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace IDTO;
using namespace DigitModified;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit-modified/data/";

int main(int argc, char* argv[]) {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);

    // Define robot model
    const std::string urdf_filename = "../Examples/Digit-modified/digit-v3-modified.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.806;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 3, 3, 3,
             3, 3, -3, 3, 3, 3, 3;
    
    // ignore all motor dynamics
    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();

    pinocchio::Data data(model);

    const double T = 0.4;
    const int N = 16;
    const TimeDiscretization time_discretization = Chebyshev;
    const int degree = 5;
    const std::string output_name = "Chebyshev-N16";

    GaitParameters gp;
    gp.swingfoot_midstep_z_des = 0.27;
    gp.swingfoot_begin_x_des = -0.25;
    gp.swingfoot_begin_y_des = 0.40;
    gp.swingfoot_end_x_des = -0.25;
    gp.swingfoot_end_y_des = -0.40;                     

    std::ifstream initial_guess(filepath + "initial-digit-modified-Bezier.txt");

    // double temp = 0;
    // std::vector<double> z_array;
    // while (initial_guess >> temp) {
    //     z_array.push_back(temp);
    // }
    // initial_guess.close();
    // Eigen::VectorXd z(z_array.size());
    // for (int i = 0; i < z_array.size(); i++) {
    //     z(i) = z_array[i];
    // }
    Eigen::VectorXd z = Eigen::VectorXd::Zero(NUM_INDEPENDENT_JOINTS * (degree + 1) + NUM_JOINTS + NUM_DEPENDENT_JOINTS);

    // // add disturbance to initial guess
    // Eigen::VectorXd disturbance(z_array.size());
    // if (argc > 1) {
    //     char* end = nullptr;
    //     std::srand((unsigned int)std::strtoul(argv[1], &end, 10));
    //     disturbance.setRandom();
    //     disturbance = disturbance * (0.2 * z.norm()) / disturbance.norm();
    // }
    // else {
    //     disturbance.setZero();
    // }

    // z = z + disturbance;
    
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

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
	app->Options()->SetNumericValue("max_wall_time", 60.0);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
    app->Options()->SetIntegerValue("max_iter", 2000);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    app->Options()->SetStringValue("nlp_scaling_method", "none");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("point_perturbation_radius", 1e-2);
    // // app->Options()->SetIntegerValue("derivative_test_first_index", 168);
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-4);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    // Run ipopt to solve the optimization problem
    double solve_time = 0;
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        // solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        solve_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        std::cout << "Total solve time: " << solve_time << " seconds.\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }
    
    // Print the solution
    // if (mynlp->solution.size() == mynlp->numVars) {
    //     std::ofstream solution(filepath + "solution-digit-modified-Bezier-" + output_name + ".txt");
    //     solution << std::setprecision(20);
    //     for (int i = 0; i < mynlp->numVars; i++) {
    //         solution << mynlp->solution[i] << std::endl;
    //     }
    //     solution.close();

    //     std::ofstream trajectory(filepath + "trajectory-digit-modified-Bezier-" + output_name + ".txt");
    //     trajectory << std::setprecision(20);
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->cidPtr_->q(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->cidPtr_->v(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->cidPtr_->a(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->cidPtr_->tau(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->cidPtr_->lambda(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     trajectory.close();
    // }

    return 0;
}
