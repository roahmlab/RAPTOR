#include "DualKinovaOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main() {
    // Define robot model
    const std::string urdf_filename1 = "../Robots/kinova-gen3/gen3_2f85_fixed.urdf";
    pinocchio::Model model1;
    pinocchio::urdf::buildModel(urdf_filename1, model1);
    model1.gravity.linear()(2) = GRAVITY;

    const std::string urdf_filename2 = "../Robots/kinova-gen3/gen3_2f85_fixed_the_other_side.urdf";
    pinocchio::Model model2;
    pinocchio::urdf::buildModel(urdf_filename2, model2);
    model2.gravity.linear()(2) = GRAVITY;

    // Define obstacles
    std::vector<Eigen::Vector3d> boxCenters = {
        Eigen::Vector3d(0.7, 0.0, 0.0), // floor
        Eigen::Vector3d(0.53, 0.49, 0.56), // back wall
        Eigen::Vector3d(-0.39, -0.82, 0.56), // bar near the control
        Eigen::Vector3d(-0.39, 0.42, 0.56), // bar between 10 and 20 change to wall
        Eigen::Vector3d(0.7, 0.0, 1.12), // ceiling
    };
    std::vector<Eigen::Vector3d> boxOrientation = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0)
    };
    std::vector<Eigen::Vector3d> boxSize = {
        Eigen::Vector3d(2.0, 2.0, 0.01),
        Eigen::Vector3d(2.0, 0.08, 1.12),
        Eigen::Vector3d(0.08, 0.08, 1.12),
        Eigen::Vector3d(0.08, 0.08, 1.12),
        Eigen::Vector3d(2.0, 2.0, 0.05)
    };

    // Define trajectories
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model1.nv + model2.nv);
    Eigen::VectorXd qT = Eigen::VectorXd::Zero(model1.nv + model2.nv);

    q0 << 0.6873863,   1.16653997,  2.11789636, -1.02120659,  0.8033111,   0.99861891, -1.39448917,
          0.6873863,   1.16653997,  2.11789636, -1.02120659,  0.8033111,   0.99861891, -1.39448917;

    qT << 0.05490901,  1.18798053,  2.03709998, -1.01553996, -1.45212933,  0.94367354, 0.46641401,
          0.05490901,  1.18798053,  2.03709998, -1.01553996, -1.45212933,  0.94367354, 0.46641401;

    const double T = 2.0;
    const int N = 100;
    const int degree = 1;

    // Define initial guess
    Eigen::VectorXd z = Eigen::VectorXd::Zero(model1.nv * degree * 3 * 2);

    Eigen::VectorXd qdiff = (qT - q0) / (degree + 1);
    for (int i = 0; i < degree; i++) {
        z.head(model1.nv * degree * 3).segment(i * model1.nv * 3, model1.nv) = q0.head(model1.nv) + qdiff.head(model1.nv) * (i + 1);
        z.tail(model1.nv * degree * 3).segment(i * model1.nv * 3, model1.nv) = q0.tail(model1.nv) + qdiff.tail(model1.nv) * (i + 1);
    }

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model1.nq);
    joint_limits_buffer.setConstant(0.0);
    Eigen::VectorXd velocity_limits_buffer(model1.nq);
    velocity_limits_buffer.setConstant(0.0);
    Eigen::VectorXd torque_limits_buffer(model1.nq);
    torque_limits_buffer.setConstant(0.0);

    // Initialize Kinova optimizer
    SmartPtr<DualKinovaOptimizer> mynlp = new DualKinovaOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              model1,
                              model2,
                              q0,
                              qT,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer,
                              true,
                              0.0,
                              0.01);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    app->Options()->SetNumericValue("obj_scaling_factor", 100.0);
	app->Options()->SetNumericValue("max_wall_time", 100.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 500);
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
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        // std::ofstream solution("solution-kinova.txt");
        // solution << std::setprecision(20);
        // for (int i = 0; i < mynlp->numVars; i++) {
        //     solution << mynlp->solution[i] << std::endl;
        // }
        // solution.close();

        std::ofstream trajectory("trajectory-dual-kinova.txt");
        trajectory << std::setprecision(20);
        for (int i = 0; i < N; i++) {
            trajectory << mynlp->kinovaOptPtr1_->trajPtr_->q(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr2_->trajPtr_->q(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr1_->trajPtr_->q_d(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr2_->trajPtr_->q_d(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr1_->trajPtr_->q_dd(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr2_->trajPtr_->q_dd(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr1_->idPtr_->tau(i).transpose() << ' ';
            trajectory << mynlp->kinovaOptPtr2_->idPtr_->tau(i).transpose() << std::endl;
        }
        trajectory.close();

        std::ofstream tapered_capsules("tapered-capsules.txt");
        const KinovaCustomizedConstraints* kccPtr1_ = dynamic_cast<
            KinovaCustomizedConstraints*>(mynlp->kinovaOptPtr1_->constraintsPtrVec_.back().get());
        const KinovaCustomizedConstraints* kccPtr2_ = dynamic_cast<
            KinovaCustomizedConstraints*>(mynlp->kinovaOptPtr2_->constraintsPtrVec_.back().get());
        const auto& tapered_capsules1 = kccPtr1_->tapered_capsules;
        const auto& tapered_capsules2 = kccPtr2_->tapered_capsules;
        for(Index i = 0; i < N; i++){
            const size_t tc1_begin_index = tapered_capsules1[2].first;
            const size_t tc1_end_index   = tapered_capsules1[2].second;
            const size_t tc2_begin_index = tapered_capsules2[2].first;
            const size_t tc2_end_index   = tapered_capsules2[2].second;

            const Eigen::Vector3d& tc1_sphere_1 = kccPtr1_->sphere_centers_copy(tc1_begin_index, i);
            const Eigen::Vector3d& tc1_sphere_2 = kccPtr1_->sphere_centers_copy(tc1_end_index, i);
            const Eigen::Vector3d& tc2_sphere_1 = kccPtr2_->sphere_centers_copy(tc2_begin_index, i);
            const Eigen::Vector3d& tc2_sphere_2 = kccPtr2_->sphere_centers_copy(tc2_end_index, i);

            tapered_capsules << tc1_sphere_1.transpose() << ' ' << tc1_sphere_2.transpose() << ' ' << tc2_sphere_1.transpose() << ' ' << tc2_sphere_2.transpose() << std::endl;
        }
        tapered_capsules.close();
    }

    return 0;
}
