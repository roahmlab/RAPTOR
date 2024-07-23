#include "KinovaWaitrOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "CustomizedInverseDynamics.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main() {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);
    
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/WaitrUnsafe/kinova_grasp.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    const int actual_model_nq = model.nq - 1;

    model.gravity.linear()(2) = -9.81;
    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3, 
             0; // the last joint is a fixed joint

    // Define obstacles
    const int num_obstacles = 2;
    std::vector<Eigen::Vector3d> boxCenters;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters[i] << 10, 0, 0.3 * i;
        boxOrientation[i] << 0, 0, 0;
        boxSize[i] << 0.1, 0.1, 0.1;
    }

    // Define trajectories
    WaitrTrajectoryParameters atp;
    atp.q0 = Eigen::VectorXd::Zero(actual_model_nq);
    atp.q_d0 = Eigen::VectorXd::Zero(actual_model_nq);
    atp.q_dd0 = Eigen::VectorXd::Zero(actual_model_nq);

    const double T = 1;
    const int N = 16;
    const int degree = WAITR_BEZIER_CURVE_DEGREE;

    // Define contact surface parameters
    contactSurfaceParams csp;
    // TODO: define contact surface parameters

    // Define target
    Eigen::VectorXd qdes(actual_model_nq);
    qdes.setConstant(1.0);
    const int tplan_n = N / 2;

    // Define initial guess
    Eigen::VectorXd z(actual_model_nq);
    z.setZero();

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(actual_model_nq);
    joint_limits_buffer.setConstant(0.0);
    Eigen::VectorXd velocity_limits_buffer(actual_model_nq);
    velocity_limits_buffer.setConstant(0.0);
    Eigen::VectorXd torque_limits_buffer(actual_model_nq);
    torque_limits_buffer.setConstant(0.0);

    // std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<WaitrBezierCurves>(T, 
    //                                                N, 
    //                                                actual_model_nq, 
    //                                                Chebyshev, 
    //                                                atp);
    // CustomizedInverseDynamics id(model, jtype, trajPtr_);

    // Eigen::VectorXd x = -5 * Eigen::VectorXd::Ones(actual_model_nq);
    // auto start = std::chrono::high_resolution_clock::now();
    // id.compute(x, true);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << "Total compute time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds.\n";
    // const Eigen::MatrixXd& pv_pz = id.ptau_pz(10);
    // Eigen::MatrixXd pv_pz_ref = Eigen::MatrixXd::Zero(pv_pz.rows(), pv_pz.cols());
    // for (int i = 0; i < actual_model_nq; i++) {
    //     Eigen::VectorXd x_plus = x;
    //     x_plus(i) += 1e-8;
    //     id.compute(x_plus, false);
    //     const Eigen::VectorXd v_plus = id.tau(10);
    //     Eigen::VectorXd x_minus = x;
    //     x_minus(i) -= 1e-8;
    //     id.compute(x_minus, false);
    //     const Eigen::VectorXd v_minus = id.tau(10);
    //     pv_pz_ref.col(i) = (v_plus - v_minus) / 2e-8;
    // }
    // // std::cout << "pv_pz: " << pv_pz << std::endl;
    // // std::cout << std::endl;
    // // std::cout << "pv_pz_ref: " << pv_pz_ref << std::endl;
    // std::cout << pv_pz - pv_pz_ref << std::endl;

    // std::cout << std::setprecision(20);
    // for (int i = 0; i < N; i++) {
    //     std::cout << "tau(" << i << ") = " << id.tau(i).transpose() << std::endl;
    // }
    // return 0;

    // Initialize Kinova optimizer
    SmartPtr<KinovaWaitrOptimizer> mynlp = new KinovaWaitrOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              model,
                              jtype,
                              atp,
                              csp,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              qdes,
                              tplan_n,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    // app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
	app->Options()->SetNumericValue("max_wall_time", 0.2);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 50);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

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
        solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Total solve time: " << solve_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        std::ofstream solution("solution-kinova-waitr.txt");
        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }
        solution.close();

        std::ofstream trajectory("trajectory-kinova-waitr.txt");
        trajectory << std::setprecision(20);
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->trajPtr_->q(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->trajPtr_->q_d(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->trajPtr_->q_dd(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->idPtr_->tau(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        trajectory.close();
    }

    return 0;
}
