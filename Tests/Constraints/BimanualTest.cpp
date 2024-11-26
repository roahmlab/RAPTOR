#define BOOST_TEST_MODULE BimanualTest
#include <boost/test/included/unit_test.hpp>

#include "ArmourBezierCurves.h"
#include "PZDynamics.h"
#include "TaperedCapsuleCollision.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(BimanualTest)

// Test collision distances for each case against the expected values

// t = 0, 0 > u > 1
BOOST_AUTO_TEST_CASE(CollisionCase1){
    Eigen::Vector3d p11(3.0, 0.0, 0.0);
    Eigen::Vector3d p12(0.0, 0.0, 0.0);
    Eigen::Vector3d p21(6.0, 3.0, 0.0);
    Eigen::Vector3d p22(6.0, -3.0, 0.0);
    double r11(1.0);
    double r12(0.5);
    double r21(1.0);
    double r22(0.5);

    TaperedCapsuleCollision<2> collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1.23957), 1e-4);
}

// t = 1, 0 > u > 1
BOOST_AUTO_TEST_CASE(CollisionCase2){
    Eigen::Vector3d p12(3.0, 0.0, 0.0);
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p21(6.0, 3.0, 0.0);
    Eigen::Vector3d p22(6.0, -3.0, 0.0);
    double r12(1.0);
    double r11(0.5);
    double r21(1.0);
    double r22(0.5);

    TaperedCapsuleCollision<2> collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1.23957), 1e-4);
}

// u = 0, 0 > t > 1
BOOST_AUTO_TEST_CASE(CollisionCase3){
    Eigen::Vector3d p21(3.0, 0.0, 0.0);
    Eigen::Vector3d p22(0.0, 0.0, 0.0);
    Eigen::Vector3d p11(6.0, 3.0, 0.0);
    Eigen::Vector3d p12(6.0, -3.0, 0.0);
    double r21(1.0);
    double r22(0.5);
    double r11(1.0);
    double r12(0.5);

    TaperedCapsuleCollision<2> collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1.23957), 1e-4);
}

// u = 1, 0 > t > 1
BOOST_AUTO_TEST_CASE(CollisionCase4){
    Eigen::Vector3d p22(3.0, 0.0, 0.0);
    Eigen::Vector3d p21(0.0, 0.0, 0.0);
    Eigen::Vector3d p11(6.0, 3.0, 0.0);
    Eigen::Vector3d p12(6.0, -3.0, 0.0);
    double r22(1.0);
    double r21(0.5);
    double r11(1.0);
    double r12(0.5);

    TaperedCapsuleCollision<2> collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1.23957), 1e-4);
}

// u = 0, t = 0
BOOST_AUTO_TEST_CASE(CollisionCase5){
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p12(-3.0, 0.0, 0.0);
    Eigen::Vector3d p21(2.0, 0.0, 0.0);
    Eigen::Vector3d p22(5.0, 0.0, 0.0);
    double r11(0.5);
    double r12(1.5);
    double r21(0.5);
    double r22(1.5);

    TaperedCapsuleCollision<7> collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1), 1e-4);
}

// u = 0, t = 1
BOOST_AUTO_TEST_CASE(CollisionCase6){
    Eigen::Vector3d p11(0.0, -3.0, 0.0);
    Eigen::Vector3d p12(0.0, 0.0, 0.0);
    Eigen::Vector3d p21(0.0, 2.0, 0.0);
    Eigen::Vector3d p22(0.0, 5.0, 0.0);
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);
    double distance = 0.0;

    TaperedCapsuleCollision<2> collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(0.5), 1e-4);
}

// u = 1, t = 0
BOOST_AUTO_TEST_CASE(CollisionCase7){
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p12(0.0, 0.0, -3.0);
    Eigen::Vector3d p21(0.0, 0.0, 5.0);
    Eigen::Vector3d p22(0.0, 0.0, 2.0);
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(0.75);
    double distance = 0.0;

    TaperedCapsuleCollision<2> collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(0.75), 1e-4);
}

// u = 1, t = 1
BOOST_AUTO_TEST_CASE(CollisionCase8){
    Eigen::Vector3d p11(0.0, 0.0, -3.0);
    Eigen::Vector3d p12(0.0, 0.0, 0.0);
    Eigen::Vector3d p21(0.0, 0.0, 5.0);
    Eigen::Vector3d p22(0.0, 0.0, 2.0);
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(0.75);
    double distance = 0.0;

    TaperedCapsuleCollision<2> collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(0.25), 1e-4);
}

// 0 < u < 1, 0 < t < 1
BOOST_AUTO_TEST_CASE(CollisionCase9){
    Eigen::Vector3d p11(5.0, 2.0, 0.0);
    Eigen::Vector3d p12(-5.0, 2.0, 0.0);
    Eigen::Vector3d p21(0.0, -2.0, 5.0);
    Eigen::Vector3d p22(0.0, -2.0, -5.0);
    double r11(1.0);
    double r12(1.5);
    double r21(1.0);
    double r22(1.5);
    double distance = 0.0;

    TaperedCapsuleCollision<2> collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1.4900), 1e-4);
}

BOOST_AUTO_TEST_CASE(GeneralGradient){
    std::srand(time(NULL));
    const double r11[5] = {0.5, 1, 1.5, 2, 5};
    const double r12[5] = {0.5, 1, 1.5, 2, 5};
    const double r21[5] = {0.5, 1, 1.5, 2, 5};
    const double r22[5] = {0.5, 1, 1.5, 2, 5};
    for(int r = 0; r<5; r++){
        for(int i = 0; i<1000; i++){
            Eigen::Vector3d p11 = Eigen::Vector3d::Random()*20;
            Eigen::Vector3d p12 = Eigen::Vector3d::Random()*20;
            Eigen::Vector3d p21 = Eigen::Vector3d::Random()*20;
            Eigen::Vector3d p22 = Eigen::Vector3d::Random()*20;
    
            Eigen::Matrix<double,3,7> p11_pz = Eigen::Matrix<double,3,7>::Random();
            Eigen::Matrix<double,3,7> p12_pz = Eigen::Matrix<double,3,7>::Random();
            Eigen::Matrix<double,3,7> p21_pz = Eigen::Matrix<double,3,7>::Random();
            Eigen::Matrix<double,3,7> p22_pz = Eigen::Matrix<double,3,7>::Random();

            TaperedCapsuleCollision<7> collider;
            Eigen::Vector<double,7> dist_grad;
            double distance = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11[r], r12[r], r21[r], r22[r], dist_grad);
            
            Eigen::Vector<double,7> z(7);
            z << 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5;
            double distance_delta = collider.computeDistance(p11+p11_pz*z, p12+p12_pz*z, p21+p21_pz*z, p22+p22_pz*z,r11[r], r12[r], r21[r], r22[r]);
            auto analytic_grad = (dist_grad.transpose()*z);
            Eigen::MatrixXd numerical_grad(1,1);
            numerical_grad << distance_delta-distance;

            BOOST_CHECK_SMALL((numerical_grad-analytic_grad).norm() , 1e-7);
        }
    }
}

BOOST_AUTO_TEST_CASE(KinovaCollisionCheck){
    #define Interval RAPTOR::Kinova::Armour::Interval
    #define Vec3 Eigen::Vector3d
    #define getCenter RAPTOR::Kinova::Armour::getCenter
        // Initialize
    // read robot model and info
    // const std::string robot_model_file = "../Robots/kinova-gen3/kinova.urdf";
    // const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithoutGripperInfo.yaml";
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova_grasp.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaSuctionCup.yaml";
    const std::shared_ptr<RAPTOR::Kinova::Armour::RobotInfo> robotInfoPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::RobotInfo>(robot_model_file, robot_info_file);

    // turn off friction for validation
    robotInfoPtr_->model.friction.setZero();

    Eigen::VectorXd q0(7);
    q0 << 0.0, 0.0, 0.0, 3.14, 0.0, 0.0, 0.0;

    Eigen::VectorXd dq0(7);
    dq0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    const double duration = 2.0;

    const Eigen::VectorXd k_center = Eigen::VectorXd::Zero(robotInfoPtr_->num_motors);
    const Eigen::VectorXd k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr_->num_motors);

    std::shared_ptr<RAPTOR::Kinova::Armour::BezierCurveInterval> trajPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::BezierCurveInterval>(
            q0, dq0, dq0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_,
            10);
    
    const int num_time_steps = trajPtr_->num_time_steps;
    std::cout << "Num time steps: " << num_time_steps << std::endl;

    std::shared_ptr<RAPTOR::Kinova::Armour::PZDynamics> dynPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::PZDynamics>(robotInfoPtr_, trajPtr_);

    std::shared_ptr<RAPTOR::TaperedCapsuleCollision<7>> tccPtr = 
        std::make_shared<RAPTOR::TaperedCapsuleCollision<7>>();

    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;

    // bimanual self collision constraints
    const int arm_1_capsule_num = robotInfoPtr_->num_capsules;
    const int arm_2_capsule_num = robotInfoPtr_->num_capsules;

    const double x_val = 1;
    const double* x;
    x = &x_val;

    std::cout << "Arm 1 Capsule Num: " << arm_1_capsule_num << std::endl;
    std::cout << "Arm 2 Capsule Num: " << arm_2_capsule_num << std::endl;

    const int check_steps = std::min(10, num_time_steps);

    for (int i = 0; i< check_steps; i++){
        for (int arm_1_index = 0; arm_1_index<arm_1_capsule_num-2; arm_1_index++){
            std::string sphere_name = robotInfoPtr_->tc_spheres[arm_1_index*2];
            pinocchio::FrameIndex frame_id = 
                robotInfoPtr_->model.getFrameId(sphere_name);

            auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();
            Interval x_res = PZsphere(0).slice(x);
            Interval y_res = PZsphere(1).slice(x);
            Interval z_res = PZsphere(2).slice(x);

            Vec3 tc1_sphere_1;
            tc1_sphere_1 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);
            // TODO: get sphere names from TC index
            sphere_name = robotInfoPtr_->tc_spheres[arm_1_index*2+1];
            frame_id = 
                robotInfoPtr_->model.getFrameId(sphere_name);

            PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();
            x_res = PZsphere(0).slice(x);
            y_res = PZsphere(1).slice(x);
            z_res = PZsphere(2).slice(x);

            Vec3 tc1_sphere_2;
            tc1_sphere_2 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);

            double tc1_sphere_1_radius = dynPtr_->sphere_radii(arm_1_index, i);
            double tc1_sphere_2_radius = dynPtr_->sphere_radii(arm_1_index+1, i);

            for (int arm_2_index = arm_1_index+2; arm_2_index<arm_2_capsule_num; arm_2_index++){
                std::string sphere_name2_1 = robotInfoPtr_->tc_spheres[arm_2_index*2];
                pinocchio::FrameIndex frame_id2 = 
                    robotInfoPtr_->model.getFrameId(sphere_name2_1);

                auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id2].translation();
                Interval x_res2 = PZsphere(0).slice(x);
                Interval y_res2 = PZsphere(1).slice(x);
                Interval z_res2 = PZsphere(2).slice(x);

                Vec3 tc2_sphere_1;
                tc2_sphere_1 << getCenter(x_res2), 
                                getCenter(y_res2), 
                                getCenter(z_res2);

                std::string sphere_name2_2 = robotInfoPtr_->tc_spheres[arm_2_index*2+1];
                frame_id2 = 
                    robotInfoPtr_->model.getFrameId(sphere_name2_2);

                PZsphere = dynPtr_->data_sparses[i].oMf[frame_id2].translation();
                x_res = PZsphere(0).slice(x);
                y_res = PZsphere(1).slice(x);
                z_res = PZsphere(2).slice(x);

                Vec3 tc2_sphere_2;
                tc2_sphere_2 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);

                const double tc2_sphere_1_radius = dynPtr_->sphere_radii(arm_2_index, i);
                const double tc2_sphere_2_radius = dynPtr_->sphere_radii(arm_2_index+1, i);
                double distance = tccPtr->computeDistance(tc1_sphere_1, tc1_sphere_2, tc2_sphere_1, tc2_sphere_2,
                                                    tc1_sphere_1_radius, tc1_sphere_2_radius, tc2_sphere_1_radius, tc2_sphere_2_radius);
                if(arm_1_index == 1 && arm_2_index == 3){
                    BOOST_CHECK(distance > 0);
                }
                else{
                    BOOST_CHECK(distance < 0);
                }
                
                std::cout << "Index1:" << arm_1_index << " Index2:" << arm_2_index << " Distance: " << distance << std::endl;
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(KinovaCollisionCheckZeroAngle){
    #define Interval RAPTOR::Kinova::Armour::Interval
    #define Vec3 Eigen::Vector3d
    #define getCenter RAPTOR::Kinova::Armour::getCenter
        // Initialize
    // read robot model and info
    // const std::string robot_model_file = "../Robots/kinova-gen3/kinova.urdf";
    // const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithoutGripperInfo.yaml";
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova_grasp.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaSuctionCup.yaml";
    const std::shared_ptr<RAPTOR::Kinova::Armour::RobotInfo> robotInfoPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::RobotInfo>(robot_model_file, robot_info_file);

    // turn off friction for validation
    robotInfoPtr_->model.friction.setZero();

    Eigen::VectorXd q0(7);
    q0 << 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0;

    Eigen::VectorXd dq0(7);
    dq0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    const double duration = 2.0;

    const Eigen::VectorXd k_center = Eigen::VectorXd::Zero(robotInfoPtr_->num_motors);
    const Eigen::VectorXd k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr_->num_motors);

    std::shared_ptr<RAPTOR::Kinova::Armour::BezierCurveInterval> trajPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::BezierCurveInterval>(
            q0, dq0, dq0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_,
            10);
    
    const int num_time_steps = trajPtr_->num_time_steps;

    std::shared_ptr<RAPTOR::Kinova::Armour::PZDynamics> dynPtr_ = 
        std::make_shared<RAPTOR::Kinova::Armour::PZDynamics>(robotInfoPtr_, trajPtr_);

    std::shared_ptr<RAPTOR::TaperedCapsuleCollision<7>> tccPtr = 
        std::make_shared<RAPTOR::TaperedCapsuleCollision<7>>();

    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();

    // bimanual self collision constraints
    const int arm_1_capsule_num = robotInfoPtr_->num_capsules;
    const int arm_2_capsule_num = robotInfoPtr_->num_capsules;

    const double x_val = 1;
    const double* x;
    x = &x_val;

    const int check_steps = std::min(10, num_time_steps);

    for (int i = 0; i< check_steps; i++){
        int num_checks = 0;
        for (int arm_1_index = 0; arm_1_index<arm_1_capsule_num-2; arm_1_index++){
            std::string sphere_name = robotInfoPtr_->tc_spheres[arm_1_index*2];
            pinocchio::FrameIndex frame_id = 
                robotInfoPtr_->model.getFrameId(sphere_name);

            auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();
            Interval x_res = PZsphere(0).slice(x);
            Interval y_res = PZsphere(1).slice(x);
            Interval z_res = PZsphere(2).slice(x);

            Vec3 tc1_sphere_1;
            tc1_sphere_1 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);
            // TODO: get sphere names from TC index
            sphere_name = robotInfoPtr_->tc_spheres[arm_1_index*2+1];
            frame_id = 
                robotInfoPtr_->model.getFrameId(sphere_name);

            PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();
            x_res = PZsphere(0).slice(x);
            y_res = PZsphere(1).slice(x);
            z_res = PZsphere(2).slice(x);

            Vec3 tc1_sphere_2;
            tc1_sphere_2 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);

            double tc1_sphere_1_radius = dynPtr_->sphere_radii(arm_1_index, i);
            double tc1_sphere_2_radius = dynPtr_->sphere_radii(arm_1_index+1, i);

            for (int arm_2_index = arm_1_index+2; arm_2_index<arm_2_capsule_num; arm_2_index++){
                std::string sphere_name2_1 = robotInfoPtr_->tc_spheres[arm_2_index*2];
                pinocchio::FrameIndex frame_id2 = 
                    robotInfoPtr_->model.getFrameId(sphere_name2_1);

                auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id2].translation();
                Interval x_res2 = PZsphere(0).slice(x);
                Interval y_res2 = PZsphere(1).slice(x);
                Interval z_res2 = PZsphere(2).slice(x);

                Vec3 tc2_sphere_1;
                tc2_sphere_1 << getCenter(x_res2), 
                                getCenter(y_res2), 
                                getCenter(z_res2);

                std::string sphere_name2_2 = robotInfoPtr_->tc_spheres[arm_2_index*2+1];
                frame_id2 = 
                    robotInfoPtr_->model.getFrameId(sphere_name2_2);

                PZsphere = dynPtr_->data_sparses[i].oMf[frame_id2].translation();
                x_res = PZsphere(0).slice(x);
                y_res = PZsphere(1).slice(x);
                z_res = PZsphere(2).slice(x);

                Vec3 tc2_sphere_2;
                tc2_sphere_2 << getCenter(x_res), 
                                getCenter(y_res), 
                                getCenter(z_res);

                const double tc2_sphere_1_radius = dynPtr_->sphere_radii(arm_2_index, i);
                const double tc2_sphere_2_radius = dynPtr_->sphere_radii(arm_2_index+1, i);
                double distance = tccPtr->computeDistance(tc1_sphere_1, tc1_sphere_2, tc2_sphere_1, tc2_sphere_2,
                                                    tc1_sphere_1_radius, tc1_sphere_2_radius, tc2_sphere_1_radius, tc2_sphere_2_radius);
                BOOST_CHECK(distance > 1e-2);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
