#define BOOST_TEST_MODULE IntervalMomentumRegressorTest
#include <boost/test/included/unit_test.hpp>

#include "IntervalMomentumRegressor.h"
#include "MomentumRegressor.h"
#include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(IntervalMomentumRegressorTest)

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
    std::shared_ptr<TrajectoryData> trajPtr_ = 
        std::make_shared<TrajectoryData>(T, N, model.nv);

    // Initialize IntervalMomentumRegressor
    IntervalMomentumRegressor regressor_id(model, trajPtr_);
    MomentumRegressor regressor_id2(model, trajPtr_);

    // This is just a placeholder, TrajectoryData does not use this
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(model.nv);

    // Compute inverse dynamics using IntervalMomentumRegressor
    regressor_id.compute(z, false);
    regressor_id2.compute(z, false);
    auto CTv = intervalDoubleMatrixMultiply(regressor_id.Y_CTv, regressor_id.phi);

    // Compute inverse dynamics using pinocchio::rnea and compare the results
    trajPtr_->compute(z, false);
    for (int i = 0; i < N; i++) {
        Eigen::VectorXd tau_pinocchio = pinocchio::rnea(
            model, data, 
            trajPtr_->q(i), 
            Eigen::VectorXd::Zero(model.nv), 
            trajPtr_->q_d(i));
        pinocchio::computeCoriolisMatrix(model, data, trajPtr_->q(i), trajPtr_->q_d(i));
        Eigen::VectorXd CTv_pinocchio = data.C.transpose() * trajPtr_->q_d(i);

        for (int j = 0; j < model.nv; j++) {
            BOOST_CHECK_SMALL(regressor_id.tau(i)(j).upper() - regressor_id.tau(i)(j).lower(), 1e-10);
            BOOST_CHECK_SMALL(regressor_id.tau(i)(j).lower() - tau_pinocchio(j), 1e-6);
            BOOST_CHECK_SMALL(CTv(i * model.nv + j).upper() - CTv(i * model.nv + j).lower(), 1e-10);
            BOOST_CHECK_SMALL(CTv(i * model.nv + j).lower() - CTv_pinocchio(j), 1e-6);
        }
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

    // Disable gravity
    model.gravity.setZero();

    // Create a trajectory
    int N = 2;  // number of time steps
    double T = 10.0;  // total time
    std::shared_ptr<TrajectoryData> trajPtr_ = 
        std::make_shared<TrajectoryData>(T, N, model.nv);

    // Initialize IntervalMomentumRegressor
    IntervalMomentumRegressor regressor_id(model, trajPtr_);

    // This is just a placeholder, TrajectoryData does not use this
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(model.nv);

    // Compute inverse dynamics using IntervalMomentumRegressor
    regressor_id.compute(z, true);

    // Compute inverse dynamics using pinocchio::rnea and compare the results
    trajPtr_->compute(z, false);
    Eigen::MatrixXd rnea_partial_dq(model.nv, model.nv);
    Eigen::MatrixXd rnea_partial_dv(model.nv, model.nv);
    Eigen::MatrixXd rnea_partial_da(model.nv, model.nv);

    for (int i = 0; i < N; i++) {
        pinocchio::computeRNEADerivatives(model, data, 
            trajPtr_->q(i), Eigen::VectorXd::Zero(model.nv), trajPtr_->q_d(i),
            rnea_partial_dq, rnea_partial_dv, rnea_partial_da);

        rnea_partial_da.triangularView<Eigen::StrictlyLower>() = rnea_partial_da.transpose().triangularView<Eigen::StrictlyLower>();

        for (int j = 0; j < model.nv; j++) {
            for (int k = 0; k < model.nv; k++) {
                BOOST_CHECK_SMALL(regressor_id.ptau_pz(i)(j, k).upper() - regressor_id.ptau_pz(i)(j, k).lower(), 1e-10);
            }

            for (int k = 0; k < model.nv; k++) {
                BOOST_CHECK_SMALL((regressor_id.ptau_pz(i)(j, k).lower() - rnea_partial_dq(j, k)), 1e-6);
                BOOST_CHECK_SMALL((regressor_id.ptau_pz(i)(j, k + model.nv).lower() - rnea_partial_da(j, k)), 1e-6);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(ComputeTestWithNoise) {
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
    SensorNoiseInfo sensor_noise(model.nv);
    sensor_noise.position_error << Utils::deg2rad(0.02),
                                   Utils::deg2rad(0.02),
                                   Utils::deg2rad(0.02),
                                   Utils::deg2rad(0.02),
                                   Utils::deg2rad(0.011),
                                   Utils::deg2rad(0.011),
                                   Utils::deg2rad(0.011); // from Kinova official support, resolution of joint encoders
    sensor_noise.velocity_error = 10 * sensor_noise.position_error;
    sensor_noise.acceleration_error = 10 * sensor_noise.position_error;   

    int N = 128;  // number of time steps
    double T = 10.0;  // total time
    std::shared_ptr<TrajectoryData> trajPtr_ = 
        std::make_shared<TrajectoryData>(T, N, model.nv, true, sensor_noise);

    // Initialize IntervalMomentumRegressor
    
    IntervalMomentumRegressor regressor_id(model, trajPtr_);

    // This is just a placeholder, TrajectoryData does not use this
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(model.nv);

    // Compute inverse dynamics using IntervalMomentumRegressor
    regressor_id.compute(z, true);
    auto CTv = intervalDoubleMatrixMultiply(regressor_id.Y_CTv, regressor_id.phi);

    // Compute inverse dynamics using pinocchio::rnea and compare the results
    trajPtr_->compute(z, false);
    for (int i = 0; i < N; i++) {
        Eigen::VectorXd tau_pinocchio = pinocchio::rnea(
            model, data, 
            trajPtr_->q(i), 
            Eigen::VectorXd::Zero(model.nv), 
            trajPtr_->q_d(i));
        pinocchio::computeCoriolisMatrix(model, data, trajPtr_->q(i), trajPtr_->q_d(i));
        Eigen::VectorXd CTv_pinocchio = data.C.transpose() * trajPtr_->q_d(i);

        for (int j = 0; j < model.nv; j++) {
            BOOST_CHECK_LE(regressor_id.tau(i)(j).lower(), tau_pinocchio(j) + 1e-14);
            BOOST_CHECK_GE(regressor_id.tau(i)(j).upper(), tau_pinocchio(j) - 1e-14);
            BOOST_CHECK_LE(CTv(i * model.nv + j).lower(), CTv_pinocchio(j) + 1e-14);
            BOOST_CHECK_GE(CTv(i * model.nv + j).upper(), CTv_pinocchio(j) - 1e-14);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()