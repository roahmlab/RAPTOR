#define BOOST_TEST_MODULE BoxCollisionAvoidanceTest
#include <boost/test/included/unit_test.hpp>

#include "BoxCollisionAvoidance.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(BoxCollisionAvoidanceTest)

Eigen::Vector3d randomVectorOnSphere(double radius) {
    // Generate random angles
    double theta = 2 * M_PI * ((double) rand() / RAND_MAX); // Azimuthal angle
    double phi = acos(2 * ((double) rand() / RAND_MAX) - 1); // Polar angle

    // Convert spherical coordinates to Cartesian coordinates
    double x = radius * sin(phi) * cos(theta);
    double y = radius * sin(phi) * sin(theta);
    double z = radius * cos(phi);

    return Eigen::Vector3d(x, y, z);
}

// test gradient
BOOST_AUTO_TEST_CASE(GradientTest)
{
    std::srand(std::time(nullptr));

    std::vector<Eigen::Vector3d> boxCenter;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenter.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    boxOrientation.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    boxSize.push_back(Eigen::Vector3d(1.0, 1.0, 1.0));

    BoxCollisionAvoidance bca(boxCenter, boxOrientation, boxSize);

    int varLength = 5;
    const double radii[3] = {1.2, 1.5, 2.0};

    for (int r = 0; r < 3; r++) {
        const double radius = radii[r];
        for (int i = 0; i < 100; i++) {
            Eigen::Vector3d point = randomVectorOnSphere(radius);
            Eigen::MatrixXd ppoint_pz = Eigen::MatrixXd::Random(3, varLength);

            bca.computeDistance(point, ppoint_pz);

            const Eigen::VectorXd J_analytical = bca.pdistances_pz.row(bca.minimumDistanceIndex);
            Eigen::VectorXd J_numerical = Eigen::VectorXd::Zero(varLength);
            for (int j = 0; j < varLength; j++) {
                Eigen::Vector3d point_plus = point + 1e-8 * ppoint_pz.col(j);
                bca.computeDistance(point_plus);
                const double f_plus = bca.minimumDistance;

                Eigen::Vector3d point_minus = point - 1e-8 * ppoint_pz.col(j);
                bca.computeDistance(point_minus);
                const double f_minus = bca.minimumDistance;

                J_numerical(j) = (f_plus - f_minus) / 2e-8;
            }

            BOOST_CHECK_SMALL((J_analytical - J_numerical).norm(), 1e-5);
        }
    }
}

// test hessian
BOOST_AUTO_TEST_CASE(HessianTest){
    std::srand(std::time(nullptr));

    std::vector<Eigen::Vector3d> boxCenter;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenter.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    boxOrientation.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    boxSize.push_back(Eigen::Vector3d(1.0, 1.0, 1.0));

    BoxCollisionAvoidance bca(boxCenter, boxOrientation, boxSize);

    int varLength = 5;
    const double radii[3] = {1.2, 1.5, 2.0};

    for (int r = 0; r < 3; r++) {
        const double radius = radii[r];
        for (int i = 0; i < 100; i++) {
            Eigen::Vector3d point = randomVectorOnSphere(radius);
            Eigen::MatrixXd ppoint_pz = Eigen::MatrixXd::Random(3, varLength);
            Eigen::Array<Eigen::MatrixXd, 3, 1> ppoint_pz_pz;
            for (int j = 0; j < 3; j++) {
                ppoint_pz_pz(j) = Eigen::MatrixXd::Random(varLength, varLength);
            }

            bca.computeDistance(point, ppoint_pz, ppoint_pz_pz);

            const Eigen::MatrixXd H_analytical = bca.pdistances_pz_pz(bca.minimumDistanceIndex);
            Eigen::MatrixXd H_numerical = Eigen::MatrixXd::Zero(varLength, varLength);
            for (int j = 0; j < varLength; j++) {
                Eigen::Vector3d point_plus = point + 1e-8 * ppoint_pz.col(j);
                Eigen::MatrixXd ppoint_pz_plus = ppoint_pz;
                for (int k = 0; k < 3; k++) {
                    ppoint_pz_plus.row(k) = ppoint_pz.row(k) + 1e-8 * ppoint_pz_pz(k).row(j);
                }
                bca.computeDistance(point_plus, ppoint_pz_plus);
                const Eigen::VectorXd J_plus = bca.pdistances_pz.row(bca.minimumDistanceIndex);
                
                Eigen::Vector3d point_minus = point - 1e-8 * ppoint_pz.col(j);
                Eigen::MatrixXd ppoint_pz_minus = ppoint_pz;
                for (int k = 0; k < 3; k++) {
                    ppoint_pz_minus.row(k) = ppoint_pz.row(k) - 1e-8 * ppoint_pz_pz(k).row(j);
                }
                bca.computeDistance(point_minus, ppoint_pz_minus);
                const Eigen::VectorXd J_minus = bca.pdistances_pz.row(bca.minimumDistanceIndex);

                const Eigen::VectorXd H_numerical_row = (J_plus - J_minus) / 2e-8;

                BOOST_CHECK_SMALL((H_analytical.row(j) - H_numerical_row.transpose()).norm(), 1e-5);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
