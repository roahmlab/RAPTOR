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

// BOOST_AUTO_TEST_CASE(GradientCase1){
//     Eigen::Vector3d p11(-1,1,1);
//     Eigen::Vector3d p12(1,-2,3);
//     Eigen::Vector3d p21(0,-2,-1);
//     Eigen::Vector3d p22(0,-3,-1);

//     Eigen::Matrix<double,3,2> p11_pz(3, 2);
//     p11_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p12_pz(3, 2);
//     p12_pz << 0,0,0,0,0.25,1;
//     Eigen::Matrix<double,3,2> p21_pz(3, 2);
//     p21_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p22_pz(3, 2);
//     p22_pz << 0,0,0,0,0.25,1;

//     double r11(0.5);
//     double r12(1.0);
//     double r21(0.5);
//     double r22(1.0);
//     double distance = 0.0;

//     TaperedCapsuleCollision collider;
//     Eigen::Vector<double,2> dist_grad;
//     distance = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);

//     BOOST_CHECK_SMALL(distance-(2.0653), 1e-4);

//     Eigen::Vector2d z(2);
//     z << 1e-4, 1e-4;
//     double distance_delta = collider.computeDistance(p11+p11_pz*z, p12+p12_pz*z, p21+p21_pz*z, p22+p22_pz*z,r11, r12, r21, r22);
//     auto analytic_grad = (dist_grad.transpose()*z);
//     Eigen::MatrixXd numerical_grad(1,1);
//     numerical_grad << distance_delta-distance;
//     BOOST_CHECK_SMALL((numerical_grad-analytic_grad).norm()*1e5 , 1e-4);
// }

// BOOST_AUTO_TEST_CASE(GradientCase2){
//     Eigen::Vector3d p11(0,-2,-1);
//     Eigen::Vector3d p12(0,-3,-1);
//     Eigen::Vector3d p21(-1,1,1);
//     Eigen::Vector3d p22(1,-2,3);

//     Eigen::Matrix<double,3,2> p11_pz(3, 2);
//     p11_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p12_pz(3, 2);
//     p12_pz << 0,0,0,0,0.25,1;
//     Eigen::Matrix<double,3,2> p21_pz(3, 2);
//     p21_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p22_pz(3, 2);
//     p22_pz << 0,0,0,0,0.25,1;

//     double r11(0.5);
//     double r12(1.0);
//     double r21(0.5);
//     double r22(1.0);
//     double distance = 0.0;

//     TaperedCapsuleCollision collider;
//     Eigen::Vector<double,2> dist_grad;
//     distance = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);

//     BOOST_CHECK_SMALL(distance-(2.0653), 1e-4);
    
//     Eigen::Vector2d z(2);
//     z << 1e-4, 1e-4;
//     double distance_delta = collider.computeDistance(p11+p11_pz*z, p12+p12_pz*z, p21+p21_pz*z, p22+p22_pz*z,r11, r12, r21, r22);
//     auto analytic_grad = (dist_grad.transpose()*z);
//     Eigen::MatrixXd numerical_grad(1,1);
//     numerical_grad << distance_delta-distance;
//     BOOST_CHECK_SMALL((numerical_grad-analytic_grad).norm()*1e5 , 1e-4);
// }

// BOOST_AUTO_TEST_CASE(GradientCase3){
//     Eigen::Vector3d p11(1,1,2);
//     Eigen::Vector3d p12(1,3,3);
//     Eigen::Vector3d p21(0,-1,-1);
//     Eigen::Vector3d p22(0,-3,-5);

//     Eigen::Matrix<double,3,2> p11_pz(3, 2);
//     p11_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p12_pz(3, 2);
//     p12_pz << 0,0,0,0,0.25,1;
//     Eigen::Matrix<double,3,2> p21_pz(3, 2);
//     p21_pz << 1,0.5,0,0,0,0;
//     Eigen::Matrix<double,3,2> p22_pz(3, 2);
//     p22_pz << 0,0,0,0,0.25,1;

//     double r11(0.5);
//     double r12(1.0);
//     double r21(0.5);
//     double r22(1.0);

//     TaperedCapsuleCollision collider;
//     Eigen::Vector<double,2> dist_grad;

//     double distance = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
    
//     Eigen::Vector2d z(2);
//     z << 1e-4, 1e-4;
//     double distance_delta = collider.computeDistance(p11+p11_pz*z, p12+p12_pz*z, p21+p21_pz*z, p22+p22_pz*z,r11, r12, r21, r22);
//     auto analytic_grad = (dist_grad.transpose()*z);
//     Eigen::MatrixXd numerical_grad(1,1);
//     numerical_grad << distance_delta-distance;
//     BOOST_CHECK_SMALL((numerical_grad-analytic_grad).norm() , 1e-4);
// }

BOOST_AUTO_TEST_SUITE_END()
