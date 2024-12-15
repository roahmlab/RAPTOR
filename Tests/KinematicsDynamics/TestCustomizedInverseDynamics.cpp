#define BOOST_TEST_MODULE  CustomizedInverseDynamicsTest
#include <boost/test/included/unit_test.hpp>

#include "pinocchio/algorithm/model.hpp"

#include "CustomizedInverseDynamics.h"
#include "BezierCurves.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(InverseDynamicsTestSuite)
 // Test without fixed joints
BOOST_AUTO_TEST_CASE(test_inverse_dynamics_without_fixed_joints)
{
    // define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova_grasp.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.armature.setZero();
    model.damping.setZero();
    model.friction.setZero();

    pinocchio::Data data(model);

    // set joint configurations
    std::shared_ptr<Trajectories> trajPtr_ = 
        std::make_shared<BezierCurves>(
            2.0, 1, model.nq, TimeDiscretization::Uniform, 5);

    std::srand(std::time(nullptr));
    Eigen::VectorXd z = Eigen::VectorXd::Random(trajPtr_->varLength);
    trajPtr_->compute(z);

    for (int i = 0; i < trajPtr_->N; i++) {
        pinocchio::rnea(model, data, 
                        trajPtr_->q(0), trajPtr_->q_d(0), trajPtr_->q_dd(0));
    }

    // compute inverse dynamics using RAPTOR
    std::shared_ptr<CustomizedInverseDynamics> cidPtr = 
        std::make_shared<CustomizedInverseDynamics>(
            model, trajPtr_);
    cidPtr->compute(z, false);

    // check the difference
    BOOST_CHECK_SMALL((data.tau - cidPtr->tau(0)).norm(), 1e-10);

    for (int i = 0; i < model.nv; i++) {
        BOOST_CHECK_SMALL((data.f[i + 1].angular() - cidPtr->f(i).head(3)).norm(), 1e-10);
        BOOST_CHECK_SMALL((data.f[i + 1].linear() - cidPtr->f(i).tail(3)).norm(), 1e-10);
    }
}

// Test with fixed joints
BOOST_AUTO_TEST_CASE(test_inverse_dynamics_with_fixed_joints)
{
    // define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova_grasp.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();
  
    std::shared_ptr<Trajectories> trajPtr_ = 
        std::make_shared<BezierCurves>(
            Eigen::VectorXd::LinSpaced(5, 0, 1), model.nq, 5);

    Eigen::VectorXi jtype = convertPinocchioJointType(model);    
    jtype(jtype.size() - 1) = 0; // fix the last joint

    pinocchio::Model model_reduced;
    std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id = {(pinocchio::JointIndex)model.nv};
    pinocchio::buildReducedModel(model, list_of_joints_to_lock_by_id, Eigen::VectorXd::Zero(model.nv), model_reduced);
    pinocchio::Data data_reduced(model_reduced);

    Eigen::VectorXd z = Eigen::VectorXd::Random(trajPtr_->varLength);

    trajPtr_ = std::make_shared<BezierCurves>(
        Eigen::VectorXd::LinSpaced(5, 0, 1), model_reduced.nq, 5);
    trajPtr_->compute(z);

    // compute inverse dynamics using pinocchio
    for (int i = 0; i < trajPtr_->N; i++) {
        pinocchio::rnea(model_reduced, data_reduced, 
                        trajPtr_->q(0).head(model_reduced.nq), 
                        trajPtr_->q_d(0).head(model_reduced.nv), 
                        trajPtr_->q_dd(0).head(model_reduced.nv));
    }
    // compute inverse dynamics using RAPTOR
    std::shared_ptr<CustomizedInverseDynamics> cidPtr = std::make_shared<CustomizedInverseDynamics>(
        model, trajPtr_, jtype);
    cidPtr->compute(z, true);

    //check the error
    BOOST_CHECK_SMALL((data_reduced.tau - cidPtr->tau(0)).norm(), 1e-10);
}
BOOST_AUTO_TEST_SUITE_END()