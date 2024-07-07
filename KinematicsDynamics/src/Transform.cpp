#include "Transform.h"

namespace IDTO {

Transform::Transform() {
    R.setIdentity();
    p.setZero();
    ifDerivative = false;
} 

Transform::Transform(const Mat3& R_in, 
                     const Vec3& p_in, 
                     const bool i_in) :
    R(R_in),
    p(p_in), 
    ifDerivative(i_in) {
}

Transform::Transform(const Vec3& p_in) :
    p(p_in) {
    R.setIdentity();
    ifDerivative = false;
}

Transform::Transform(const int jtype, 
                     const double theta, 
                     const int order) {
    if (order == 0) { // Transform of rotation matrix
        R = Mat3::Identity();
    }
    else if (order == 1) { // First order derivative of transform of rotation matrix
        R = Mat3::Zero();
    }
    else if (order == 2) { // Second order derivative of transform of rotation matrix
        R = Mat3::Zero();
    }
    else if (order == 3) { // Third order derivative of transform of rotation matrix
        R = Mat3::Zero();
    }
    else {
        throw std::runtime_error("spatial.cpp: Transform(): Unsupported differentiate order!");
    }
    p = Vec3::Zero(3); // translation vector
    
    ifDerivative = (order > 0);

    if (order == 0) { // original transformation
        if (jtype == 1) {
            R << 1, 0, 0,
                 0, cos(theta), -sin(theta),
                 0, sin(theta), cos(theta);
        }
        else if (jtype == -1) {
            R << 1, 0, 0,
                 0, cos(theta), sin(theta),
                 0, -sin(theta), cos(theta);
        }
        else if (jtype == 2) {
            R << cos(theta), 0, sin(theta),
                 0, 1, 0,
                 -sin(theta), 0, cos(theta);
        }
        else if (jtype == -2) {
            R << cos(theta), 0, -sin(theta),
                 0, 1, 0,
                 sin(theta), 0, cos(theta);
        }
        else if (jtype == 3) {
            R << cos(theta), -sin(theta), 0,
                 sin(theta), cos(theta), 0,
                 0, 0, 1;
        }
        else if (jtype == -3) {
            R << cos(theta), sin(theta), 0,
                 -sin(theta), cos(theta), 0,
                 0, 0, 1;
        }
        else if (jtype == 4) {
            p << theta, 0, 0;
        }
        else if (jtype == 5) {  
            p << 0, theta, 0;
        }
        else if (jtype == 6) {
            p << 0, 0, theta;
        }
        else if (jtype == -4) {
            p << -theta, 0, 0;
        }
        else if (jtype == -5) {  
            p << 0, -theta, 0;
        }
        else if (jtype == -6) {
            p << 0, 0, -theta;
        }
        else {
            throw std::runtime_error("spatial.cpp: Transform(): Unrecognized jtype!");
        }
    }
    else if (order == 1) { // first order derivative
        if (jtype == 1) {
            R << 0, 0, 0,
                 0, -sin(theta), -cos(theta),
                 0, cos(theta), -sin(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, -sin(theta), cos(theta),
                 0, -cos(theta), -sin(theta);
        }
        else if (jtype == 2) {
            R << -sin(theta), 0, cos(theta),
                 0, 0, 0,
                 -cos(theta), 0, -sin(theta);
        }
        else if (jtype == -2) {
            R << -sin(theta), 0, -cos(theta),
                 0, 0, 0,
                 cos(theta), 0, -sin(theta);
        }
        else if (jtype == 3) {
            R << -sin(theta), -cos(theta), 0,
                 cos(theta), -sin(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << -sin(theta), cos(theta), 0,
                 -cos(theta), -sin(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == 4) {
            p << 1, 0, 0;
        }
        else if (jtype == 5) {
            p << 0, 1, 0;
        }
        else if (jtype == 6) {
            p << 0, 0, 1;
        }
        else if (jtype == -4) {
            p << -1, 0, 0;
        }
        else if (jtype == -5) {
            p << 0, -1, 0;
        }
        else if (jtype == -6) {
            p << 0, 0, -1;
        }
        else {
            throw std::runtime_error("spatial.cpp: Transform(): Unrecognized jtype!");
        }
    }
    else if (order == 2) { // Second order derivative
        if (jtype == 1) {
            R << 0, 0, 0,
                 0, -cos(theta), sin(theta),
                 0, -sin(theta), -cos(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, -cos(theta), -sin(theta),
                 0, sin(theta), -cos(theta);
        }
        else if (jtype == 2) {
            R << -cos(theta), 0, -sin(theta),
                 0, 0, 0,
                 sin(theta), 0, -cos(theta);
        }
        else if (jtype == -2) {
            R << -cos(theta), 0, sin(theta),
                 0, 0, 0,
                 -sin(theta), 0, -cos(theta);
        }
        else if (jtype == 3) {
            R << -cos(theta), sin(theta), 0,
                 -sin(theta), -cos(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << -cos(theta), -sin(theta), 0,
                 sin(theta), -cos(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == 4) {
            p << 0, 0, 0;
        }
        else if (jtype == 5) {
            p << 0, 0, 0;
        }
        else if (jtype == 6) {
            p << 0, 0, 0;
        }
        else if (jtype == -4) {
            p << 0, 0, 0;
        }
        else if (jtype == -5) {
            p << 0, 0, 0;
        }
        else if (jtype == -6) {
            p << 0, 0, 0;
        }
        else {
            throw std::runtime_error("spatial.cpp: Transform(): Unrecognized jtype!");
        }
    }
    else if (order == 3) { // Third order derivative
        if (jtype == 1) {
            R << 0, 0, 0,
                 0, sin(theta), cos(theta),
                 0, -cos(theta), sin(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, sin(theta), -cos(theta),
                 0, cos(theta), sin(theta);
        }
        else if (jtype == 2) {
            R << sin(theta), 0, -cos(theta),
                 0, 0, 0,
                 cos(theta), 0, sin(theta);
        }
        else if (jtype == -2) {
            R << sin(theta), 0, cos(theta),
                 0, 0, 0,
                 -cos(theta), 0, sin(theta);
        }
        else if (jtype == 3) {
            R << sin(theta), cos(theta), 0,
                 -cos(theta), sin(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << sin(theta), -cos(theta), 0,
                 cos(theta), sin(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == 4) {
            p << 0, 0, 0;
        }
        else if (jtype == 5) {
            p << 0, 0, 0;
        }
        else if (jtype == 6) {
            p << 0, 0, 0;
        }
        else if (jtype == -4) {
            p << 0, 0, 0;
        }
        else if (jtype == -5) {
            p << 0, 0, 0;
        }
        else if (jtype == -6) {
            p << 0, 0, 0;
        }
        else {
            throw std::runtime_error("spatial.cpp: Transform(): Unrecognized jtype!");
        }
    }
    else {
        throw std::runtime_error("spatial.cpp: Transform(): Unsupported differentiate order!");
    }
}

Transform Transform::operator=(const Transform& x) {
    R = x.R;
    p = x.p;
    ifDerivative = x.ifDerivative;
    return *this;
}

Transform Transform::operator*(const Transform& x) const {
    bool i = ifDerivative | x.ifDerivative;
    if (x.ifDerivative) return Transform(R * x.R, R * x.p, i);
    return Transform(R * x.R, R * x.p + p, i);
}

Transform Transform::operator*=(const Transform& x) {
    if (x.ifDerivative) {
        p = R * x.p;
    }
    else {
        p = R * x.p + p;
    }
    R = R * x.R;
    ifDerivative |= x.ifDerivative;
    return *this;
}

Transform Transform::operator*(const SE3& x) const {
    return Transform(R * x.rotation(), 
                     R * x.translation() + p, 
                     ifDerivative);
}

Transform Transform::operator*=(const SE3& x) {
    p = R * x.translation() + p;
    R = R * x.rotation();
    return *this;
}

Transform operator*(const pinocchio::SE3Tpl<double>& x, const Transform& sRp) {
    return Transform(x.rotation() * sRp.R, 
                     x.rotation() * sRp.p + x.translation(), 
                     sRp.ifDerivative);
}

Transform Transform::inverse() const {
    assert(!ifDerivative);
    return Transform(R.transpose(), -R.transpose() * p);
}

Eigen::Vector3d Transform::getRPY() const {
    // roll pitch yaw conversion does not make sense 
    // for the derivative of a transformation matrix
    assert(!ifDerivative);

    Eigen::Vector3d rpy;
    rpy << atan2(-R(1,2), R(2,2)),                   // roll
           HigherOrderDerivatives::safeasin(R(0,2)), // pitch
           atan2(-R(0,1), R(0,0));                   // yaw

    return rpy;
}

Eigen::Matrix<double, 6, 1> Transform::getXYZRPY() const {
    Eigen::Matrix<double, 6, 1> xyzrpy;
    xyzrpy << p, getRPY();
    return xyzrpy;
}

std::ostream& operator<<(std::ostream& os, const Transform& sRp) {
    os << "rotation:\n"    << sRp.R << std::endl 
       << "translation:\n" << sRp.p << std::endl;
    return os;
}

}  // namespace IDTO
