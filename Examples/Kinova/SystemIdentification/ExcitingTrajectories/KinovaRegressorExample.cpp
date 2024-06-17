#include "ConditionNumberOptimizer.h"

using namespace IDTO;
using namespace Kinova;
using namespace Ipopt;

int main(int argc, char* argv[]) {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);

    const std::string regroupMatrixFileName = "../Examples/Kinova/SystemIdentification/RegroupMatrix.csv";
    
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.81;
    model.friction.setZero();
    // model.damping.setZero();
    // model.rotorInertia.setZero();

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    // Define trajectory parameters
    double T = 5.0;
    int N = 100;
    int degree = 5;

    // Define initial guess
    Eigen::VectorXd z = 1.0 * Eigen::VectorXd::Random((2 * degree + 3) * model.nv);
    z.block((2 * degree + 1) * model.nv, 0, model.nv, 1) = 
        2.23 * Eigen::VectorXd::Random(model.nv);
    z.block((2 * degree + 1) * model.nv + model.nv, 0, model.nv, 1)
        = Eigen::VectorXd::Zero(model.nv);

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model.nq);
    joint_limits_buffer.setConstant(0.0);
    Eigen::VectorXd velocity_limits_buffer(model.nq);
    velocity_limits_buffer.setConstant(0.0);
    Eigen::VectorXd torque_limits_buffer(model.nq);
    torque_limits_buffer.setConstant(0.0);

    // Initialize Kinova optimizer
    SmartPtr<ConditionNumberOptimizer> mynlp = new ConditionNumberOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              model,
                              jtype,
                              regroupMatrixFileName,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
	app->Options()->SetNumericValue("max_wall_time", 40.0);
	app->Options()->SetIntegerValue("print_level", 5);
    // app->Options()->SetIntegerValue("max_iter", 50);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
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
        solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Total solve time: " << solve_time << " milliseconds.\n";
    }
    catch (int errorCode) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Re-evaluate the solution at a higher resolution and print the results.
    std::shared_ptr<Trajectories> traj = std::make_shared<FixedFrequencyFourierCurves>(T, 
                                                                                       2000, 
                                                                                       model.nv, 
                                                                                       TimeDiscretization::Uniform, 
                                                                                       degree);
    traj->compute(mynlp->solution, false);

    if (argc > 1) {
        std::ofstream position("../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/exciting-position-" + std::string(argv[1]) + ".csv");
        std::ofstream velocity("../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/exciting-velocity-" + std::string(argv[1]) + ".csv");
        std::ofstream acceleration("../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/exciting-acceleration-" + std::string(argv[1]) + ".csv");

        for (int i = 0; i < traj->N; i++) {
            position << traj->q(i).transpose() << std::endl;
            velocity << traj->q_d(i).transpose() << std::endl;
            acceleration << traj->q_dd(i).transpose() << std::endl;
        }
    }
    else {
        std::ofstream position("exciting-position.csv");
        std::ofstream velocity("exciting-velocity.csv");
        std::ofstream acceleration("exciting-acceleration.csv");

        for (int i = 0; i < traj->N; i++) {
            position << traj->q(i).transpose() << std::endl;
            velocity << traj->q_d(i).transpose() << std::endl;
            acceleration << traj->q_dd(i).transpose() << std::endl;
        }
    }


    // std::shared_ptr<Trajectories> traj = std::make_shared<FixedFrequencyFourierCurves>(1.0, 
    //                                                                                    5, 
    //                                                                                    model.nv, 
    //                                                                                    TimeDiscretization::Uniform, 
    //                                                                                    4);

    // RegressorInverseDynamics rid(model, jtype, traj);
    // InverseDynamics id(model, traj);

    // Eigen::VectorXd z = Eigen::VectorXd::Random(traj->varLength);

    // auto start = std::chrono::high_resolution_clock::now();
    // rid.compute(z, true);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Regressor Inverse Dynamics: " << duration.count() << " microseconds" << std::endl;

    // start = std::chrono::high_resolution_clock::now();
    // id.compute(z, false);
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Inverse Dynamics: " << duration.count() << " microseconds" << std::endl;

    // for (int i = 0; i < rid.tau.size(); i++) {
    //     std::cout << rid.tau(i).transpose() - id.tau(i).transpose() << std::endl;
    // }

    // // compute numerical gradient by central finite difference
    // for (int j = 0; j < traj->N; j++) {
    //     for (int i = 0; i < z.size(); i++) {
    //         Eigen::VectorXd z_plus = z;
    //         Eigen::VectorXd z_minus = z;
    //         z_plus(i) += 1e-8;
    //         z_minus(i) -= 1e-8;

    //         rid.compute(z_plus, false);
    //         Eigen::VectorXd tau_plus = rid.tau(j);

    //         rid.compute(z_minus, false);
    //         Eigen::VectorXd tau_minus = rid.tau(j);

    //         Eigen::VectorXd grad = (tau_plus - tau_minus) / 2e-8;
    //         std::cout << (grad - rid.ptau_pz(j).col(i)).norm() << std::endl;
    //         // std::cout <<  << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}