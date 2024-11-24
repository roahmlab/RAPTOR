#define BOOST_TEST_MODULE BimanualTest
#include <boost/test/included/unit_test.hpp>

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

            BOOST_CHECK_SMALL((numerical_grad-analytic_grad).norm() , 1e-8);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
