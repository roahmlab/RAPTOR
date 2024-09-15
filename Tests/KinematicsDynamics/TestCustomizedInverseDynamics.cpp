#define BOOST_TEST_MODULE  CustomizedInverseDynamicsTest
#include <boost/test/included/unit_test.hpp>

#include "pinocchio/algorithm/model.hpp"

#include "CustomizedInverseDynamics.h"
#include "BezierCurves.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(InverseDynamicsTestSuite)

BOOST_AUTO_TEST_CASE(test_inverse_dynamics_without_fixed_joints)
{
    // define robot model
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        perror("getcwd() error");
    }
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova_grasp.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();

    pinocchio::Data data(model);

    // set joint configurations
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<BezierCurves>(
            Eigen::VectorXd::LinSpaced(5, 0, 1), model.nq, 5);

    std::srand(std::time(nullptr));
    Eigen::VectorXd z = Eigen::VectorXd::Random(trajPtr->varLength);
    trajPtr->compute(z);

// Test without fixed joints
    // compute inverse dynamics using pinocchio
    // auto start_clock = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < trajPtr->N; i++) {
        pinocchio::rnea(model, data, 
                        trajPtr->q(0), trajPtr->q_d(0), trajPtr->q_dd(0));
    }
    // auto stop_clock = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    // std::cout << "Pinocchio ID: " << duration.count() << " nanoseconds" << std::endl;

    // compute inverse dynamics using RAPTOR
    std::shared_ptr<CustomizedInverseDynamics> cidPtr = 
        std::make_shared<CustomizedInverseDynamics>(
            model, trajPtr);

    // start_clock = std::chrono::high_resolution_clock::now();
    cidPtr->compute(z, false);
    // stop_clock = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    // std::cout << "RAPTOR ID: " << duration.count() << " nanoseconds" << std::endl;

    // compare the results
    // std::cout << "Pinocchio: " << data.tau.transpose() << std::endl;
    // std::cout << "RAPTOR: " << cidPtr->tau(0).transpose() << std::endl;
    BOOST_CHECK_SMALL((data.tau - cidPtr->tau(0)).norm(), 1e-10);

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
  
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<BezierCurves>(
            Eigen::VectorXd::LinSpaced(5, 0, 1), model.nq, 5);



    Eigen::VectorXi jtype = convertPinocchioJointType(model);    
    jtype(jtype.size() - 1) = 0; // fix the last joint

    pinocchio::Model model_reduced;
    std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id = {(pinocchio::JointIndex)model.nv};
    pinocchio::buildReducedModel(model, list_of_joints_to_lock_by_id, Eigen::VectorXd::Zero(model.nv), model_reduced);
    pinocchio::Data data_reduced(model_reduced);

    Eigen::VectorXd z = Eigen::VectorXd::Random(trajPtr->varLength);

    trajPtr = std::make_shared<BezierCurves>(
        Eigen::VectorXd::LinSpaced(5, 0, 1), model_reduced.nq, 5);
    trajPtr->compute(z);

    // compute inverse dynamics using pinocchio
    // start_clock = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < trajPtr->N; i++) {
        pinocchio::rnea(model_reduced, data_reduced, 
                        trajPtr->q(0).head(model_reduced.nq), 
                        trajPtr->q_d(0).head(model_reduced.nv), 
                        trajPtr->q_dd(0).head(model_reduced.nv));
    }
    // compute inverse dynamics using RAPTOR
    std::shared_ptr<CustomizedInverseDynamics> cidPtr = std::make_shared<CustomizedInverseDynamics>(
        model, trajPtr, jtype);
    cidPtr->compute(z, true);

    
    // stop_clock = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    // std::cout << "Pinocchio ID (fixed joints): " << duration.count() << " nanoseconds" << std::endl;

    // compute inverse dynamics using RAPTOR
    // std::shared_ptr<CustomizedInverseDynamics> cidPtr = std::make_shared<CustomizedInverseDynamics>(
    //     model, trajPtr, jtype);
    // // auto start = std::chrono::high_resolution_clock::now();
    // cidPtr->compute(z, true);
    // auto end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // std::cout << "RAPTOR ID (fixed joints): " << duration.count() << " nanoseconds" << std::endl;

    // // compare the results
    // std::cout << "Pinocchio (fixed joints): " << data_reduced.tau.transpose() << std::endl;
    // std::cout << "RAPTOR (fixed joints): " << cidPtr->tau(0).transpose() << std::endl;


    BOOST_CHECK_SMALL((data_reduced.tau - cidPtr->tau(0)).norm(), 1e-10);

 
}
BOOST_AUTO_TEST_SUITE_END()