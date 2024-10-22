#include "DigitSystemIdentification.h"

using namespace RAPTOR;
using namespace DigitWholeBodySysID;

int main() {
    // #ifdef NUM_THREADS
    //     Eigen::setNbThreads(NUM_THREADS);
    // #else
    //     throw std::runtime_error("macro NUM_THREADS is not defined!");
    // #endif

    // define robot model
    const char stanceLeg = 'L';

    // const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // Add left foot contact frame
    pinocchio::SE3 left_foot_placement;
    left_foot_placement.rotation().matrix() <<
        0, 1, 0,
        -0.5, 0, sin(M_PI/3),
        sin(M_PI/3), 0, 0.5;
    left_foot_placement.translation() <<
        0, -0.05456, -0.0315;
    model.addFrame(
        pinocchio::Frame(
            "left_foot",
            model.getJointId("left_toe_roll"),
            0,
            left_foot_placement,
            pinocchio::OP_FRAME
        )
    );

    // Add right foot contact frame
    pinocchio::SE3 right_foot_placement;
    right_foot_placement.rotation().matrix() <<
        0, -1, 0,
        0.5, 0, -sin(M_PI/3),
        sin(M_PI/3), 0, 0.5;
    right_foot_placement.translation() <<
        0, 0.05456, -0.0315;
    model.addFrame(
        pinocchio::Frame(
            "right_foot",
            model.getJointId("right_toe_roll"),
            0,
            right_foot_placement,
            pinocchio::OP_FRAME
        )
    );

    model.gravity.linear()(2) = -9.806;
    
    // ignore friction for now
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

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

    // read trajectory from data
    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>();

    *posDataPtr_ = Utils::initializeEigenMatrixFromFile("../data_stand_realworld/Digit_data_stand_realworld_position.txt");
    *torqueDataPtr_ = Utils::initializeEigenMatrixFromFile("../data_stand_realworld/Digit_data_stand_realworld_torque.txt");

    std::cout << "Number of samples: " << posDataPtr_->cols() << std::endl;

    // Initialize system identification optimizer
    SmartPtr<DigitSystemIdentification> mynlp = new DigitSystemIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              torqueDataPtr_);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-3);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-3);
    mynlp->constr_viol_tol = 1e-4;
	app->Options()->SetNumericValue("max_wall_time", 1000.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 500);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
    
    app->Options()->SetStringValue("jac_c_constant", "yes");

    // // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

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

    // Eigen::VectorXd residual = mynlp->A * mynlp->solution - mynlp->b;
    // for (int i = 0; i < mynlp->N; i++) {
    //     std::cout << residual.segment(i * NUM_INDEPENDENT_JOINTS, NUM_INDEPENDENT_JOINTS).transpose() << std::endl;
    // }

    const Eigen::VectorXd& lambdas = mynlp->solution.tail(mynlp->N * (NUM_DEPENDENT_JOINTS + 6));
    for (int i = 0; i < mynlp->N; i++) {
        const auto& lambda = lambdas.segment(i * (NUM_DEPENDENT_JOINTS + 6), NUM_DEPENDENT_JOINTS + 6);
        std::cout << lambda.tail(12).transpose() << std::endl;
    }

    pinocchio::Model new_model = model;
    
    for (int i = 0; i < mynlp->nontrivialLinkIds.size(); i++) {
        std::cout << model.names[mynlp->nontrivialLinkIds[i] + 1] << std::endl;
        std::cout << mynlp->solution.segment(i * 10, 10).transpose() << std::endl;
        std::cout << "New\n";
        std::cout << pinocchio::Inertia::FromDynamicParameters(
            mynlp->solution.segment(i * 10, 10)) << std::endl << std::endl;
        std::cout << "Old\n";
        std::cout << pinocchio::Inertia::FromDynamicParameters(
            mynlp->x0.segment(i * 10, 10)) << std::endl << std::endl;
        std::cout << std::endl;

        new_model.inertias[mynlp->nontrivialLinkIds[i] + 1] = 
            pinocchio::Inertia::FromDynamicParameters(mynlp->solution.segment(i * 10, 10));

        const double mass = mynlp->solution(10 * i);
        const Eigen::Vector3d& com = mynlp->solution.segment(10 * i + 1, 3);
        const pinocchio::Symmetric3Tpl<double> inertia(mynlp->solution.segment(10 * i + 4, 6));
        
        // construct LMI matrix
        Eigen::Matrix4d LMI = Eigen::Matrix4d::Zero();

        LMI.topLeftCorner<3, 3>() = 
            0.5 * inertia.matrix().trace() * Eigen::Matrix3d::Identity() - 
            inertia.matrix();

        LMI.topRightCorner<3, 1>() = com;
        LMI.bottomLeftCorner<1, 3>() = com.transpose();

        LMI(3, 3) = mass;

        std::cout << LMI << std::endl << std::endl;
    }

    pinocchio::Data new_data(new_model);
    const auto ddcPtr_ = mynlp->ddcPtr_;
    std::ofstream tau_indep("tau_indep_realworld.txt");
    std::ofstream tau_dep("tau_dep_realworld.txt");

    for (int i = 0; i < mynlp->N; i++) {
        const auto& lambda = lambdas.segment(i * (NUM_DEPENDENT_JOINTS + 6), NUM_DEPENDENT_JOINTS + 6);
        pinocchio::rnea(
            new_model, new_data, 
            posDataPtr_->col(i), 
            Eigen::VectorXd::Zero(new_model.nv), 
            Eigen::VectorXd::Zero(new_model.nv));
        ddcPtr_->get_J(posDataPtr_->col(i));
        Eigen::VectorXd tau_estimated = new_data.tau - ddcPtr_->J.transpose() * lambda;
        
        tau_indep << ddcPtr_->get_independent_vector(tau_estimated).transpose() << std::endl;
        tau_dep << ddcPtr_->get_dependent_vector(tau_estimated).transpose() << std::endl;
    }   

    tau_indep.close();
    tau_dep.close();
    
    return 0;
}