#define BOOST_TEST_MODULE RegressorInverseDynamicsTest
#include <boost/test/included/unit_test.hpp>

#include "RegressorInverseDynamics.h"
#include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "Polynomials.h"

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(RegressorInverseDynamicsTest)

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

    // Create a trajectory
    int N = 128;  // number of time steps
    double T = 10.0;  // total time
    int degree = 5;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr_ = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize RegressorInverseDynamics
    RegressorInverseDynamics regressor_id(model, trajPtr_, false);

    // Generate random joint p, v, and a
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr_->varLength);

    // Compute inverse dynamics using RegressorInverseDynamics
    regressor_id.compute(z, false);

    // Compute inverse dynamics using pinocchio::rnea and compare the results
    trajPtr_->compute(z, false);
    for (int i = 0; i < N; i++) {
        Eigen::VectorXd tau_pinocchio = pinocchio::rnea(
            model, data, 
            trajPtr_->q(i), 
            trajPtr_->q_d(i),
            trajPtr_->q_dd(i));

        BOOST_CHECK_SMALL((regressor_id.tau(i) - tau_pinocchio).norm(), 1e-6);
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
    int N = 10;  // number of time steps
    double T = 1.0;  // total time
    int degree = 3;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr_ = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize RegressorInverseDynamics
    RegressorInverseDynamics regressor_id(model, trajPtr_);

    // Generate random joint p, v, and a
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = Eigen::VectorXd::Random(trajPtr_->varLength);

    // Compute inverse dynamics using RegressorInverseDynamics
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
        BOOST_CHECK_SMALL((pY_pz(i) - J_Y_numerical).cwiseAbs().maxCoeff(), 1e-5);
    }
}

BOOST_AUTO_TEST_SUITE_END()