#define BOOST_TEST_MODULE RegressorInverseDynamicsTest
#include <boost/test/included/unit_test.hpp>
#include "RegressorInverseDynamics.h"
// #include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "Polynomials.h"

using VecX = Eigen::VectorXd;
using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(RegressorInverseDynamicsSuite)
BOOST_AUTO_TEST_CASE(RegressorInverseDynamicsAccuracy)
{
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // Disable rotor inertia, friction, and damping
    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero();

    // Create a trajectory
    int N = 5;  // number of time steps
    double T = 10.0;  // total time
    int degree = 5;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize RegressorInverseDynamics
    RegressorInverseDynamics regressor_id(model, trajPtr);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    VecX z = 2 * M_PI * VecX::Random(trajPtr->varLength).array() - M_PI;

    // Compute inverse dynamics using RegressorInverseDynamics
    regressor_id.compute(z, false);

    // Compute inverse dynamics using pinocchio::rnea
    trajPtr->compute(z, false);
    VecX tau_pinocchio = pinocchio::rnea(model, data, trajPtr->q(0), trajPtr->q_d(0), trajPtr->q_dd(0));

    // compare the results
    BOOST_CHECK_SMALL((regressor_id.tau(0) -tau_pinocchio).norm(), 1e-10);
}
BOOST_AUTO_TEST_SUITE_END()