#define BOOST_TEST_MODULE ForwardKinematicsTest
#include <boost/test/included/unit_test.hpp>
#include "ForwardKinematics.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(ForwardKinematicsSuite)
BOOST_AUTO_TEST_CASE(TestForwardKinematicsAccuracy)
{
    // Define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    ForwardKinematicsSolver fkSolver(&model);

    // set joint angles
    std::srand(std::time(nullptr));
    Eigen::VectorXd q = M_2_PI * Eigen::VectorXd::Random(model.nq);

    // compute forward kinematics using pinocchio
    pinocchio::forwardKinematics(model, data, q);

    // set the start and end joint
    int start = 0;
    int end = model.getJointId("left_toe_B");
    fkSolver.compute(start, end, q);;

    // compare the results
    Eigen::Vector3d pinocchio_translation = data.oMi[model.getJointId("left_toe_B")].translation();
    Eigen::Vector3d raptor_translation = fkSolver.getTranslation();

    // check the error
    BOOST_CHECK_SMALL((pinocchio_translation - raptor_translation).norm(), 1e-10);
}
BOOST_AUTO_TEST_SUITE_END()