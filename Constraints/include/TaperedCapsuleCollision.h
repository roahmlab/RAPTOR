#ifndef TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
#define TAPERED_CAPSULE_COLLISION_AVOIDANCE_H

#include <Eigen/Dense>
#include <iostream>
#include "Utils.h"
#include <vector>

namespace RAPTOR {

#define COLLISION_THRESHOLD 1e-4

double solve_quadratic(float a, float b, float c, int sign);

class TaperedCapsuleCollision {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TaperedCapsuleCollision() = default;

    // Destructor
    ~TaperedCapsuleCollision() = default;

    // class methods:

    double computeDistance(const Eigen::Vector3d& tc1_point_1, const Eigen::Vector3d& tc1_point_2, 
                        const Eigen::Vector3d& tc2_point_1, const Eigen::Vector3d& tc2_point_2, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2);
                         
};

}; // namespace RAPTOR

#endif // TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
