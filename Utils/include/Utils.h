#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <Eigen/Dense>

namespace IDTO {

inline double deg2rad(const double deg) {
    return deg * M_PI / 180.0;
}

inline double safeasin(const double a) {
    if (a >= 1.0) {
        return M_PI / 2.0;
    } 
    else if (a <= -1.0) {
        return -M_PI / 2.0;
    } 
    else {
        return asin(a);
    }
}

inline double wrapToPi(const double angle) {
    double res = angle;
    while (res > M_PI) {
        res -= 2.0 * M_PI;
    }
    while (res < -M_PI) {
        res += 2.0 * M_PI;
    }
    return res;
}

inline Eigen::VectorXd wrapToPi(const Eigen::VectorXd& angles) {
    Eigen::VectorXd res = angles;
    for (int i = 0; i < res.size(); i++) {
        res(i) = wrapToPi(res(i));
    }
    return res;
}

template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

inline bool ifTwoVectorEqual(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double tol = 1e-10) {
    if (a.size() != b.size()) {
        return false;
    }
    for (int i = 0; i < a.size(); i++) {
        if (fabs(a(i) - b(i)) > tol) {
            return false;
        }
    }
    return true;
}

inline Eigen::MatrixXd reshape(const Eigen::VectorXd& vec, int rows, int cols) {
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), rows, cols);
}

inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d res;
    res << 0,    -v(2), v(1),
           v(2),  0,   -v(0),
          -v(1),  v(0), 0;
    return res;
}

inline Eigen::Vector3d skew(const Eigen::Matrix3d& m) {
    Eigen::Vector3d res;
    res << m(2,1) - m(1,2), 
           m(0,2) - m(2,0), 
           m(1,0) - m(0,1);
    return res;
}

inline Eigen::Matrix<double, 6, 6> plux(const Eigen::Matrix3d& R, const Eigen::Vector3d& p) {
    Eigen::Matrix<double, 6, 6> res;
    res << R,            Eigen::MatrixXd::Zero(3, 3),
           -R * skew(p), R;
    return res;
}

}; // namespace IDTO

#endif // UTILS_H
