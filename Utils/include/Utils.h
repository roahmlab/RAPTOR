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

}; // namespace IDTO

#endif // UTILS_H
