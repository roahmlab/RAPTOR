#include "DigitSystemIdentification.h"

using namespace RAPTOR;
using namespace DigitWholeBodySysID;

int main() {
    // define robot model
    const char stanceLeg = 'L';

    // const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // model.gravity.linear()(2) = GRAVITY;
    
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
    // const std::string trajectory_filename = "../Examples/Digit/data/trajectory-digit-Bezier.txt";
    // const Eigen::MatrixXd trajectory = Utils::initializeEigenMatrixFromFile(trajectory_filename);

    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>();

    for (int data_id = 0; data_id < 16; data_id++) {
        const auto posData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_position.txt");
        const auto velData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_velocity.txt");
        const auto accData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_acceleration.txt");
        const auto torqueData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_torque.txt");

        if (posData.cols() != velData.cols() || posData.cols() != accData.cols() || posData.cols() != torqueData.cols()) {
            throw std::runtime_error("Data size mismatch!");
        }

        if (data_id == 0) {
            *posDataPtr_ = posData;
            *velDataPtr_ = velData;
            *accDataPtr_ = accData;
            *torqueDataPtr_ = torqueData;
        }
        else {
            posDataPtr_->conservativeResize(Eigen::NoChange, posDataPtr_->cols() + posData.cols());
            velDataPtr_->conservativeResize(Eigen::NoChange, velDataPtr_->cols() + velData.cols());
            accDataPtr_->conservativeResize(Eigen::NoChange, accDataPtr_->cols() + accData.cols());
            torqueDataPtr_->conservativeResize(Eigen::NoChange, torqueDataPtr_->cols() + torqueData.cols());

            posDataPtr_->block(0, posDataPtr_->cols() - posData.cols(), posData.rows(), posData.cols()) = posData;
            velDataPtr_->block(0, velDataPtr_->cols() - velData.cols(), velData.rows(), velData.cols()) = velData;
            accDataPtr_->block(0, accDataPtr_->cols() - accData.cols(), accData.rows(), accData.cols()) = accData;
            torqueDataPtr_->block(0, torqueDataPtr_->cols() - torqueData.cols(), torqueData.rows(), torqueData.cols()) = torqueData;
        }
    }

    std::cout << "Number of samples: " << posDataPtr_->cols() << std::endl;

    std::ofstream file_pos_train("../data_realworld/left_stance_segment_position_train.txt");
    file_pos_train << *posDataPtr_ << std::endl;
    file_pos_train.close();
    std::ofstream file_vel_train("../data_realworld/left_stance_segment_velocity_train.txt");
    file_vel_train << *velDataPtr_ << std::endl;
    file_vel_train.close();
    std::ofstream file_acc_train("../data_realworld/left_stance_segment_acceleration_train.txt");
    file_acc_train << *accDataPtr_ << std::endl;
    file_acc_train.close();
    std::ofstream file_tau_train("../data_realworld/left_stance_segment_torque_train.txt");
    file_tau_train << *torqueDataPtr_ << std::endl;
    file_tau_train.close();

    // Initialize system identification optimizer
    SmartPtr<DigitSystemIdentification> mynlp = new DigitSystemIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_,
                              stanceLeg);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-5);
    mynlp->constr_viol_tol = 1e-5;
	app->Options()->SetNumericValue("max_wall_time", 200.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 200);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma86");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }

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

    // load test data
    posDataPtr_->resize(0, 0);
    velDataPtr_->resize(0, 0);
    accDataPtr_->resize(0, 0);
    torqueDataPtr_->resize(0, 0);

    for (int data_id = 16; data_id < 21; data_id++) {
        const auto posData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_position.txt");
        const auto velData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_velocity.txt");
        const auto accData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_acceleration.txt");
        const auto torqueData = Utils::initializeEigenMatrixFromFile("../data_realworld/left_stance_segment_" + std::to_string(data_id) + "_torque.txt");

        if (posData.cols() != velData.cols() || posData.cols() != accData.cols() || posData.cols() != torqueData.cols()) {
            throw std::runtime_error("Data size mismatch!");
        }

        if (data_id == 16) {
            *posDataPtr_ = posData;
            *velDataPtr_ = velData;
            *accDataPtr_ = accData;
            *torqueDataPtr_ = torqueData;
        }
        else {
            posDataPtr_->conservativeResize(Eigen::NoChange, posDataPtr_->cols() + posData.cols());
            velDataPtr_->conservativeResize(Eigen::NoChange, velDataPtr_->cols() + velData.cols());
            accDataPtr_->conservativeResize(Eigen::NoChange, accDataPtr_->cols() + accData.cols());
            torqueDataPtr_->conservativeResize(Eigen::NoChange, torqueDataPtr_->cols() + torqueData.cols());

            posDataPtr_->block(0, posDataPtr_->cols() - posData.cols(), posData.rows(), posData.cols()) = posData;
            velDataPtr_->block(0, velDataPtr_->cols() - velData.cols(), velData.rows(), velData.cols()) = velData;
            accDataPtr_->block(0, accDataPtr_->cols() - accData.cols(), accData.rows(), accData.cols()) = accData;
            torqueDataPtr_->block(0, torqueDataPtr_->cols() - torqueData.cols(), torqueData.rows(), torqueData.cols()) = torqueData;
        }
    }

    SmartPtr<DigitSystemIdentification> testnlp = new DigitSystemIdentification();
    try {
	    testnlp->set_parameters(model,
                                posDataPtr_,
                                velDataPtr_,
                                accDataPtr_,
                                torqueDataPtr_,
                                stanceLeg);

        Index n, m, nnz_jac_g, nnz_h_lag;
        TNLP::IndexStyleEnum index_style;
        testnlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

        Number x[testnlp->numVars];
        Number f;

        for (int i = 0; i < testnlp->numVars; i++) {
            x[i] = mynlp->solution(i);
        }

        testnlp->eval_f(n, x, true, f);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    const auto& tau_estimated = testnlp->tau_estimated;
    Eigen::MatrixXd tau_estimated_mat = Eigen::Map<const Eigen::MatrixXd>(tau_estimated.data(), testnlp->Nact, testnlp->N);
    std::ofstream file_tau_estimated("../data_realworld/left_stance_segment_torque_estimated.txt");
    file_tau_estimated << tau_estimated_mat << std::endl; 
    file_tau_estimated.close();

    std::ofstream file_pos_test("../data_realworld/left_stance_segment_position_test.txt");
    file_pos_test << *posDataPtr_ << std::endl;
    file_pos_test.close();
    std::ofstream file_vel_test("../data_realworld/left_stance_segment_velocity_test.txt");
    file_vel_test << *velDataPtr_ << std::endl;
    file_vel_test.close();
    std::ofstream file_acc_test("../data_realworld/left_stance_segment_acceleration_test.txt");
    file_acc_test << *accDataPtr_ << std::endl;
    file_acc_test.close();
    std::ofstream file_tau_test("../data_realworld/left_stance_segment_torque_test.txt");
    file_tau_test << *torqueDataPtr_ << std::endl; 
    file_tau_test.close();

    std::cout << mynlp->solution << std::endl;

    return 0;
}