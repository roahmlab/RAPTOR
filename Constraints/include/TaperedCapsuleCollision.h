#ifndef TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
#define TAPERED_CAPSULE_COLLISION_AVOIDANCE_H

#include <Eigen/Dense>
#include <iostream>
#include "Utils.h"
#include <vector>

namespace RAPTOR {

#define COLLISION_THRESHOLD 1e-4

inline void solve_quadratic(const double a, 
                            const double b, 
                            const double c,
                            double* sol1,
                            double* sol2);

template<int factors>
class TaperedCapsuleCollision {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecF = Eigen::Vector<double, factors>;
    using Mat3F = Eigen::Matrix<double, 3, factors>;

    inline Eigen::Vector<double,factors> batchDot(const Vec3 vector, 
                                                  const Mat3F matrix);


    inline double distanceInternal(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                                   const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                                   const Mat3F& ptc1_point_1_pz, const Mat3F& ptc1_point_2_pz, 
                                   const Mat3F& ptc2_point_1_pz, const Mat3F& ptc2_point_2_pz, 
                                   const double tc1_radius_1, const double tc1_radius_2, 
                                   const double tc2_radius_1, const double tc2_radius_2,
                                   VecF& pdist_pz);

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
                           const Mat3F& ptc1_point_1_pz, const Mat3F& ptc1_point_2_pz, 
                           const Mat3F& ptc2_point_1_pz, const Mat3F& ptc2_point_2_pz, 
                           const double tc1_radius_1, const double tc1_radius_2, 
                           const double tc2_radius_1, const double tc2_radius_2,
                           VecF& pdiff_pz);
};

}; // namespace RAPTOR

#endif // TAPERED_CAPSULE_COLLISION_AVOIDANCE_H
