#include "DigitSystemIdentification.h"

using namespace RAPTOR;
using namespace Digit;

int main() {
    // define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;
    
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
    const std::string trajectory_filename = "../Examples/Digit/data/trajectory-digit-Bezier.txt";
    const Eigen::MatrixXd trajectory = Utils::initializeEigenMatrixFromFile(trajectory_filename);

    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>();

    const int N = trajectory.cols();
    *posDataPtr_ = trajectory.topRows(NUM_JOINTS);
    *velDataPtr_ = trajectory.middleRows(NUM_JOINTS, NUM_JOINTS);
    *accDataPtr_ = trajectory.middleRows(2 * NUM_JOINTS, NUM_JOINTS);
    *torqueDataPtr_ = trajectory.middleRows(3 * NUM_JOINTS, NUM_INDEPENDENT_JOINTS);

    // Initialize system identification optimizer
    SmartPtr<DigitSystemIdentification> mynlp = new DigitSystemIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_);
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
    app->Options()->SetIntegerValue("max_iter", 1000);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma57");
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

    return 0;
}