#define BOOST_TEST_MODULE MomentumRegressorTest
#include <boost/test/included/unit_test.hpp>

#include "MomentumRegressor.h"
#include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "Polynomials.h"

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(MomentumRegressorTest)

BOOST_AUTO_TEST_CASE(ComputeTest) {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // Disable rotor inertia, friction, and damping
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();
    
    // Disable gravity
    model.gravity.setZero();

    // Create a trajectory
    int N = 128;  // number of time steps
    double T = 10.0;  // total time
    int degree = 5;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize MomentumRegressor
    MomentumRegressor regressor_id(model, trajPtr);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr->varLength);

    // Compute momentums using RegressorInverseDynamics
    regressor_id.compute(z, false);

    // Compute momentums using pinocchio::rnea and compare the results
    trajPtr->compute(z, false);
    for (int i = 0; i < N; i++) {
        Eigen::VectorXd momentum_pinocchio = pinocchio::rnea(
            model, data, 
            trajPtr->q(i), 
            Eigen::VectorXd::Zero(model.nv), 
            trajPtr->q_d(i));

        BOOST_CHECK_SMALL((regressor_id.tau(i) - momentum_pinocchio).norm(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(GradientTest) {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // Disable rotor inertia, friction, and damping
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

    // Create a trajectory
    int N = 2;  // number of time steps
    double T = 10.0;  // total time
    int degree = 2;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize MomentumRegressor
    MomentumRegressor regressor_id(model, trajPtr);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr->varLength);

    // Compute momentums using RegressorInverseDynamics
    regressor_id.compute(z, true);
    const auto pY_pz = regressor_id.pY_pz;

    // Compute numerical gradient
    for (int i = 0; i < z.size(); i++) {
        Eigen::VectorXd z_plus = z;
        z_plus(i) += 1e-8;
        regressor_id.compute(z_plus, false);
        const Eigen::VectorXd tau_plus = regressor_id.tau(0);
        const Eigen::MatrixXd Y_plus = regressor_id.Y;

        Eigen::VectorXd z_minus = z;
        z_minus(i) -= 1e-8;
        regressor_id.compute(z_minus, false);
        const Eigen::VectorXd tau_minus = regressor_id.tau(0);
        const Eigen::MatrixXd Y_minus = regressor_id.Y;

        Eigen::VectorXd J_tau_numerical = (tau_plus - tau_minus) / 2e-8;
        Eigen::MatrixXd J_Y_numerical = (Y_plus - Y_minus) / 2e-8;

        BOOST_CHECK_SMALL((regressor_id.ptau_pz(0).col(i) - J_tau_numerical).norm(), 1e-5);
        BOOST_CHECK_SMALL((pY_pz(i) - J_Y_numerical).cwiseAbs().maxCoeff(), 1e-3);
    }
}

BOOST_AUTO_TEST_SUITE_END()