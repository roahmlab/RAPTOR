#include "Transform.h"

namespace RAPTOR {

Eigen::VectorXi convertPinocchioJointType(const pinocchio::ModelTpl<float>& model) {
    Eigen::VectorXi jtype(model.nq);

    for (int i = 0; i < model.nq; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint
        const std::string shortname = model.joints[pinocchio_joint_id].shortname();

        if (shortname.find('R') != std::string::npos) {
            if (shortname.find('X') != std::string::npos) {
                jtype(i) = 1;
            }
            else if (shortname.find('Y') != std::string::npos) {
                jtype(i) = 2;
            }
            else if (shortname.find('Z') != std::string::npos) {
                jtype(i) = 3;
            }
            else if (shortname.find('U') != std::string::npos) {
                // This is specific to the digit robot
                // There are 4 joints that have "<axis xyz="0 0 -1"/>"
                // But they can not be identified by urdf parser so we manually set them to be of type -3
                std::cerr << "Warning: joint " << i << " is set to be of type -3!" << std::endl;
                std::cerr << "         We assume that this would only happen for Digit-v3, but not other robots!" << std::endl;
                jtype(i) = -3;
            }
            else {
                throw std::invalid_argument("convertPinocchioJointType: invalid joint type!");
            }
        }
        else if (shortname.find('P') != std::string::npos) {
            if (shortname.find('X') != std::string::npos) {
                jtype(i) = 4;
            }
            else if (shortname.find('Y') != std::string::npos) {
                jtype(i) = 5;
            }
            else if (shortname.find('Z') != std::string::npos) {
                jtype(i) = 6;
            }
            else {
                throw std::invalid_argument("convertPinocchioJointType: invalid joint type!");
            }
        }
        else {
            throw std::invalid_argument("convertPinocchioJointType: invalid joint type!");
        }
    }

    return jtype;
}

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

Transform::Transform(const Vec3& rpy_in,
                     const Vec3& p_in) :
    p(p_in) {
    R = Eigen::AngleAxisf(rpy_in(0), Eigen::Vector3f::UnitX()) * 
        Eigen::AngleAxisf(rpy_in(1), Eigen::Vector3f::UnitY()) * 
        Eigen::AngleAxisf(rpy_in(2), Eigen::Vector3f::UnitZ());
    ifDerivative = false;
}

Transform::Transform(const Vec3& p_in) :
    p(p_in) {
    R.setIdentity();
    ifDerivative = false;
}

Transform::Transform(const int jtype, 
                     const float theta, 
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
                 0, cosf(theta), -sinf(theta),
                 0, sinf(theta), cosf(theta);
        }
        else if (jtype == -1) {
            R << 1, 0, 0,
                 0, cosf(theta), sinf(theta),
                 0, -sinf(theta), cosf(theta);
        }
        else if (jtype == 2) {
            R << cosf(theta), 0, sinf(theta),
                 0, 1, 0,
                 -sinf(theta), 0, cosf(theta);
        }
        else if (jtype == -2) {
            R << cosf(theta), 0, -sinf(theta),
                 0, 1, 0,
                 sinf(theta), 0, cosf(theta);
        }
        else if (jtype == 3) {
            R << cosf(theta), -sinf(theta), 0,
                 sinf(theta), cosf(theta), 0,
                 0, 0, 1;
        }
        else if (jtype == -3) {
            R << cosf(theta), sinf(theta), 0,
                 -sinf(theta), cosf(theta), 0,
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
                 0, -sinf(theta), -cosf(theta),
                 0, cosf(theta), -sinf(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, -sinf(theta), cosf(theta),
                 0, -cosf(theta), -sinf(theta);
        }
        else if (jtype == 2) {
            R << -sinf(theta), 0, cosf(theta),
                 0, 0, 0,
                 -cosf(theta), 0, -sinf(theta);
        }
        else if (jtype == -2) {
            R << -sinf(theta), 0, -cosf(theta),
                 0, 0, 0,
                 cosf(theta), 0, -sinf(theta);
        }
        else if (jtype == 3) {
            R << -sinf(theta), -cosf(theta), 0,
                 cosf(theta), -sinf(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << -sinf(theta), cosf(theta), 0,
                 -cosf(theta), -sinf(theta), 0,
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
                 0, -cosf(theta), sinf(theta),
                 0, -sinf(theta), -cosf(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, -cosf(theta), -sinf(theta),
                 0, sinf(theta), -cosf(theta);
        }
        else if (jtype == 2) {
            R << -cosf(theta), 0, -sinf(theta),
                 0, 0, 0,
                 sinf(theta), 0, -cosf(theta);
        }
        else if (jtype == -2) {
            R << -cosf(theta), 0, sinf(theta),
                 0, 0, 0,
                 -sinf(theta), 0, -cosf(theta);
        }
        else if (jtype == 3) {
            R << -cosf(theta), sinf(theta), 0,
                 -sinf(theta), -cosf(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << -cosf(theta), -sinf(theta), 0,
                 sinf(theta), -cosf(theta), 0,
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
                 0, sinf(theta), cosf(theta),
                 0, -cosf(theta), sinf(theta);
        }
        else if (jtype == -1) {
            R << 0, 0, 0,
                 0, sinf(theta), -cosf(theta),
                 0, cosf(theta), sinf(theta);
        }
        else if (jtype == 2) {
            R << sinf(theta), 0, -cosf(theta),
                 0, 0, 0,
                 cosf(theta), 0, sinf(theta);
        }
        else if (jtype == -2) {
            R << sinf(theta), 0, cosf(theta),
                 0, 0, 0,
                 -cosf(theta), 0, sinf(theta);
        }
        else if (jtype == 3) {
            R << sinf(theta), cosf(theta), 0,
                 -cosf(theta), sinf(theta), 0,
                 0, 0, 0;
        }
        else if (jtype == -3) {
            R << sinf(theta), -cosf(theta), 0,
                 cosf(theta), sinf(theta), 0,
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

Transform operator*(const pinocchio::SE3Tpl<float>& x, const Transform& sRp) {
    return Transform(x.rotation() * sRp.R, 
                     x.rotation() * sRp.p + x.translation(), 
                     sRp.ifDerivative);
}

Transform Transform::inverse() const {
    assert(!ifDerivative);
    return Transform(R.transpose(), -R.transpose() * p, false);
}

Eigen::Vector3f Transform::getRPY() const {
    // roll pitch yaw conversion does not make sense 
    // for the derivative of a transformation matrix
    assert(!ifDerivative);

    Eigen::Vector3f rpy;
    rpy << atan2(-R(1,2), R(2,2)),                   // roll
           HigherOrderDerivatives::safeasin(R(0,2)), // pitch
           atan2(-R(0,1), R(0,0));                   // yaw

    return rpy;
}

Eigen::Matrix<float, 6, 1> Transform::getXYZRPY() const {
    Eigen::Matrix<float, 6, 1> xyzrpy;
    xyzrpy << p, getRPY();
    return xyzrpy;
}

std::ostream& operator<<(std::ostream& os, const Transform& sRp) {
    os << "rotation:\n"    << sRp.R << std::endl 
       << "translation:\n" << sRp.p << std::endl;
    return os;
}

}  // namespace RAPTOR
