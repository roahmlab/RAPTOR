#ifndef TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
#define TAPERED_CAPSULE_COLLISION_AVOIDANCE_H

#include <Eigen/Dense>
#include <iostream>
#include "Utils.h"
#include <vector>

namespace RAPTOR {

#define COLLISION_THRESHOLD 1e-4
#define NUM_FACTORS 3

double solve_quadratic(float a, float b, float c, int sign);

using Vec3 = Eigen::Vector3d;
using MatX = Eigen::MatrixXd;
MatX batchDot(const Vec3 vector, const MatX);

double distanceInternal(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2,
                        MatX& pdist_pz);

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
    
    double computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2);

    double computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2,
                        MatX& pdiff_pz);
                         
};

}; // namespace RAPTOR

#endif // TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
