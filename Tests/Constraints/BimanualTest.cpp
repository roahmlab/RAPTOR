#define BOOST_TEST_MODULE BimanualTest
#include <boost/test/included/unit_test.hpp>

#include "TaperedCapsuleCollision.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(BimanualTest)

// test case 1

BOOST_AUTO_TEST_CASE(CollisionCase1){
    Eigen::Vector3d p11(0.0, -2.0, 3.0);
    Eigen::Vector3d p12(0.0, -2.0, 1.0);
    Eigen::Vector3d p21(0.0, -2.0, -1.0);
    Eigen::Vector3d p22(0.0, -2.0, -3.0);
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);

    TaperedCapsuleCollision collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(0.5), 1e-4);
}

BOOST_AUTO_TEST_CASE(CollisionCase2){
    Eigen::Vector3d p11(0.0, -2.0, 1.0);
    Eigen::Vector3d p12(0.0, -2.0, 3.0);
    Eigen::Vector3d p21(0.0, -2.0, -1.0);
    Eigen::Vector3d p22(0.0, -2.0, -3.0);
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);

    TaperedCapsuleCollision collider;
    double distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(1), 1e-4);
}

BOOST_AUTO_TEST_CASE(CollisionCase3){
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p12(1.0, -1.0, 0.0);
    Eigen::Vector3d p21(1.0, 2.0, 0.0);
    Eigen::Vector3d p22(2.0, 3.0, 0.0);
    double r11(1.5);
    double r12(1.0);
    double r21(1.0);
    double r22(1.5);
    double distance = 0.0;

    TaperedCapsuleCollision collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(-0.2639), 1e-4);
}

BOOST_AUTO_TEST_CASE(CollisionCase4){
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p12(1.0, -1.0, 0.0);
    Eigen::Vector3d p21(1.0, 2.0, 0.0);
    Eigen::Vector3d p22(2.0, 3.0, 0.0);
    double r11(1.0);
    double r12(1.5);
    double r21(1.5);
    double r22(1.0);
    double distance = 0.0;

    TaperedCapsuleCollision collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(-0.2657), 1e-4);
}

BOOST_AUTO_TEST_CASE(CollisionCase5){
    Eigen::Vector3d p11(0.0, 0.0, 0.0);
    Eigen::Vector3d p12(1.0, -1.0, 0.0);
    Eigen::Vector3d p21(1.0, 2.0, 0.0);
    Eigen::Vector3d p22(2.0, 3.0, 0.0);
    double r11(1.0);
    double r12(1.5);
    double r21(1.5);
    double r22(1.0);
    double distance = 0.0;

    TaperedCapsuleCollision collider;
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    BOOST_CHECK_SMALL(distance-(-0.2657), 1e-4);
}

BOOST_AUTO_TEST_CASE(GradientCase1){
    Eigen::Vector3d p11(-1,1,1);
    Eigen::Vector3d p12(1,-2,3);
    Eigen::Vector3d p21(0,-2,-1);
    Eigen::Vector3d p22(0,-3,-1);

    Eigen::MatrixXd p11_pz(3, 2);
    p11_pz << 1,0.5,0,0,0,0;
    Eigen::MatrixXd p12_pz(3, 2);
    p12_pz << 1,0.5,0,0,0,0;
    Eigen::MatrixXd p21_pz(3, 2);
    p21_pz << 0,0,0,0,0.25,1;
    Eigen::MatrixXd p22_pz(3, 2);
    p22_pz << 0,0,0,0,0.25,1;

    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);
    double distance = 0.0;

    TaperedCapsuleCollision collider;
    Eigen::MatrixXd dist_grad;
    distance = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);

    std::cout << "Distance: " << distance << "\n";
    std::cout << "Dist Grad: " << dist_grad << "\n";

    BOOST_CHECK_SMALL(distance-(2.0653), 1e-4);

    Eigen::MatrixXd correctGrad(1,2);
    correctGrad << 0.0337, 0.2425;
    BOOST_CHECK_SMALL((dist_grad-correctGrad).norm(), 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
