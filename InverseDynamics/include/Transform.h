#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace IDTO {

// 6x6 spatial transform matrix
class Transform {
public:
    using SE3 = pinocchio::SE3Tpl<double>;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;

    // Constructor
    Transform();

    // Constructor
    Transform(const Mat3& R_in, const Vec3& p_in, const bool i_in = false) : 
        R(R_in) , p(p_in), ifDerivative(i_in) {}

    // Constructor
    Transform(const Vec3& p_in) : 
        p(p_in) {
        R.setIdentity();
    }

    // Constructor
    Transform(const int jtype, const double theta, const int order = 0);

    // Destructor
    ~Transform() = default;

    // class methods:
    Transform operator=(const Transform& x);

    Transform operator*(const Transform& x) const;

    Transform operator*=(const Transform& x);

    Transform operator*(const SE3& x) const;

    Transform operator*=(const SE3& x);

    Transform inverse() const;

    // class members:
    Mat3 R; // rotation matrix
    Vec3 p; // translation vector

    // this is like the right bottom entry of a 4x4 transformation matrix.
    // usually it should be 1, but for derivatives, it should be 0.
    bool ifDerivative = false; 
};

std::ostream& operator<<(std::ostream& os, const Transform& sRp);

}  // namespace IDTO

#endif