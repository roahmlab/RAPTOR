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

}; // namespace IDTO

#endif // UTILS_H
