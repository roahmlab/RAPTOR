#include "Transform.h"

namespace IDTO {

Transform::Transform() {
    std::memset(data, 0, sizeof(data));
    data[0] = 1.0;
    data[4] = 1.0;
    data[8] = 1.0;
    ifDerivative = false;
} 

Transform::Transform(const Mat3& R_in, 
                     const Vec3& p_in, 
                     const bool i_in) : 
    ifDerivative(i_in) {
    data[0] = R_in(0, 0);
    data[1] = R_in(0, 1);
    data[2] = R_in(0, 2);
    data[3] = R_in(1, 0);
    data[4] = R_in(1, 1);
    data[5] = R_in(1, 2);
    data[6] = R_in(2, 0);
    data[7] = R_in(2, 1);
    data[8] = R_in(2, 2);
    data[9] = p_in(0);
    data[10] = p_in(1);
    data[11] = p_in(2);
}

Transform::Transform(const Vec3& p_in) {
    std::memset(data, 0, sizeof(data));
    data[0] = 1.0;
    data[4] = 1.0;
    data[8] = 1.0;
    data[9] = p_in(0);
    data[10] = p_in(1);
    data[11] = p_in(2);
}

Transform::Transform(const int jtype, 
                     const double theta, 
                     const int order) {
    Mat3 R;
    Vec3 p;
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

    data[0] = R(0, 0);
    data[1] = R(0, 1);
    data[2] = R(0, 2);
    data[3] = R(1, 0);
    data[4] = R(1, 1);
    data[5] = R(1, 2);
    data[6] = R(2, 0);
    data[7] = R(2, 1);
    data[8] = R(2, 2);
    data[9] = p(0);
    data[10] = p(1);
    data[11] = p(2);
}

Transform Transform::operator=(const Transform& x) {
    std::memcpy(data, x.data, sizeof(data));
    ifDerivative = x.ifDerivative;
    return *this;
}

Transform Transform::operator*(const Transform& x) const {
    Transform result;
    result.ifDerivative = ifDerivative | x.ifDerivative;

    if (x.ifDerivative) {
        result.data[0] = data[0] * x.data[0] + data[1] * x.data[3] + data[2] * x.data[6];
        result.data[1] = data[0] * x.data[1] + data[1] * x.data[4] + data[2] * x.data[7];
        result.data[2] = data[0] * x.data[2] + data[1] * x.data[5] + data[2] * x.data[8];
        result.data[3] = data[3] * x.data[0] + data[4] * x.data[3] + data[5] * x.data[6];
        result.data[4] = data[3] * x.data[1] + data[4] * x.data[4] + data[5] * x.data[7];
        result.data[5] = data[3] * x.data[2] + data[4] * x.data[5] + data[5] * x.data[8];
        result.data[6] = data[6] * x.data[0] + data[7] * x.data[3] + data[8] * x.data[6];
        result.data[7] = data[6] * x.data[1] + data[7] * x.data[4] + data[8] * x.data[7];
        result.data[8] = data[6] * x.data[2] + data[7] * x.data[5] + data[8] * x.data[8];
        result.data[9] = data[0] * x.data[9] + data[1] * x.data[10] + data[2] * x.data[11];
        result.data[10] = data[3] * x.data[9] + data[4] * x.data[10] + data[5] * x.data[11];
        result.data[11] = data[6] * x.data[9] + data[7] * x.data[10] + data[8] * x.data[11];
    }
    else {
        result.data[0] = data[0] * x.data[0] + data[1] * x.data[3] + data[2] * x.data[6];
        result.data[1] = data[0] * x.data[1] + data[1] * x.data[4] + data[2] * x.data[7];
        result.data[2] = data[0] * x.data[2] + data[1] * x.data[5] + data[2] * x.data[8];
        result.data[3] = data[3] * x.data[0] + data[4] * x.data[3] + data[5] * x.data[6];
        result.data[4] = data[3] * x.data[1] + data[4] * x.data[4] + data[5] * x.data[7];
        result.data[5] = data[3] * x.data[2] + data[4] * x.data[5] + data[5] * x.data[8];
        result.data[6] = data[6] * x.data[0] + data[7] * x.data[3] + data[8] * x.data[6];
        result.data[7] = data[6] * x.data[1] + data[7] * x.data[4] + data[8] * x.data[7];
        result.data[8] = data[6] * x.data[2] + data[7] * x.data[5] + data[8] * x.data[8];
        result.data[9] = data[0] * x.data[9] + data[1] * x.data[10] + data[2] * x.data[11] + data[9];
        result.data[10] = data[3] * x.data[9] + data[4] * x.data[10] + data[5] * x.data[11] + data[10];
        result.data[11] = data[6] * x.data[9] + data[7] * x.data[10] + data[8] * x.data[11] + data[11];
    }

    return result;
}

Transform Transform::operator*=(const Transform& x) {
    *this = *this * x;
    return *this;
}

Transform Transform::operator*(const SE3& x) const {
    const Mat3& xR = x.rotation();
    const Vec3& xp = x.translation();

    Transform result;
    result.ifDerivative = ifDerivative;

    result.data[0] = data[0] * xR(0, 0) + data[1] * xR(1, 0) + data[2] * xR(2, 0);
    result.data[1] = data[0] * xR(0, 1) + data[1] * xR(1, 1) + data[2] * xR(2, 1);
    result.data[2] = data[0] * xR(0, 2) + data[1] * xR(1, 2) + data[2] * xR(2, 2);
    result.data[3] = data[3] * xR(0, 0) + data[4] * xR(1, 0) + data[5] * xR(2, 0);
    result.data[4] = data[3] * xR(0, 1) + data[4] * xR(1, 1) + data[5] * xR(2, 1);
    result.data[5] = data[3] * xR(0, 2) + data[4] * xR(1, 2) + data[5] * xR(2, 2);
    result.data[6] = data[6] * xR(0, 0) + data[7] * xR(1, 0) + data[8] * xR(2, 0);
    result.data[7] = data[6] * xR(0, 1) + data[7] * xR(1, 1) + data[8] * xR(2, 1);
    result.data[8] = data[6] * xR(0, 2) + data[7] * xR(1, 2) + data[8] * xR(2, 2);
    result.data[9] = data[0] * xp(0) + data[1] * xp(1) + data[2] * xp(2) + data[9];
    result.data[10] = data[3] * xp(0) + data[4] * xp(1) + data[5] * xp(2) + data[10];
    result.data[11] = data[6] * xp(0) + data[7] * xp(1) + data[8] * xp(2) + data[11];

    return result;
}

Transform operator*(const pinocchio::SE3Tpl<double>& x, const Transform& sRp) {
    const Eigen::Matrix3d& xR = x.rotation();
    const Eigen::Vector3d& xp = x.translation();

    Transform result;
    result.ifDerivative = sRp.ifDerivative;

    result.data[0] = xR(0, 0) * sRp.data[0] + xR(0, 1) * sRp.data[3] + xR(0, 2) * sRp.data[6];
    result.data[1] = xR(0, 0) * sRp.data[1] + xR(0, 1) * sRp.data[4] + xR(0, 2) * sRp.data[7];
    result.data[2] = xR(0, 0) * sRp.data[2] + xR(0, 1) * sRp.data[5] + xR(0, 2) * sRp.data[8];
    result.data[3] = xR(1, 0) * sRp.data[0] + xR(1, 1) * sRp.data[3] + xR(1, 2) * sRp.data[6];
    result.data[4] = xR(1, 0) * sRp.data[1] + xR(1, 1) * sRp.data[4] + xR(1, 2) * sRp.data[7];
    result.data[5] = xR(1, 0) * sRp.data[2] + xR(1, 1) * sRp.data[5] + xR(1, 2) * sRp.data[8];
    result.data[6] = xR(2, 0) * sRp.data[0] + xR(2, 1) * sRp.data[3] + xR(2, 2) * sRp.data[6];
    result.data[7] = xR(2, 0) * sRp.data[1] + xR(2, 1) * sRp.data[4] + xR(2, 2) * sRp.data[7];
    result.data[8] = xR(2, 0) * sRp.data[2] + xR(2, 1) * sRp.data[5] + xR(2, 2) * sRp.data[8];
    result.data[9] = xR(0, 0) * sRp.data[9] + xR(0, 1) * sRp.data[10] + xR(0, 2) * sRp.data[11] + xp(0);
    result.data[10] = xR(1, 0) * sRp.data[9] + xR(1, 1) * sRp.data[10] + xR(1, 2) * sRp.data[11] + xp(1);
    result.data[11] = xR(2, 0) * sRp.data[9] + xR(2, 1) * sRp.data[10] + xR(2, 2) * sRp.data[11] + xp(2);

    return result;
}

Transform Transform::operator*=(const SE3& x) {
    *this = *this * x;
    return *this;
}

Transform Transform::inverse() const {
    assert(!ifDerivative);
   
    Transform result;
    result.ifDerivative = false;

    result.data[0] = data[0];
    result.data[1] = data[3];
    result.data[2] = data[6];
    result.data[3] = data[1];
    result.data[4] = data[4];
    result.data[5] = data[7];
    result.data[6] = data[2];
    result.data[7] = data[5];
    result.data[8] = data[8];
    result.data[9] = -(result.data[0] * data[9] + result.data[1] * data[10] + result.data[2] * data[11]);
    result.data[10] = -(result.data[3] * data[9] + result.data[4] * data[10] + result.data[5] * data[11]);
    result.data[11] = -(result.data[6] * data[9] + result.data[7] * data[10] + result.data[8] * data[11]);

    return result;
}

std::ostream& operator<<(std::ostream& os, const Transform& sRp) {
    os << "R = [" << sRp.data[0] << ", " << sRp.data[1] << ", " << sRp.data[2] << "; " 
                  << sRp.data[3] << ", " << sRp.data[4] << ", " << sRp.data[5] << "; " 
                  << sRp.data[6] << ", " << sRp.data[7] << ", " << sRp.data[8] << "], "
       << "p = [" << sRp.data[9] << ", " << sRp.data[10] << ", " << sRp.data[11] << "]\n";
    
    return os;
}

Eigen::Vector3d Transform::getTranslation() const {
    return Eigen::Vector3d(data[9], data[10], data[11]);
}

}  // namespace IDTO
