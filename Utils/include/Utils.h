#ifndef UTILS_H
#define UTILS_H

#include <cmath>

namespace IDTO {

inline double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

inline double safeasin(double a) {
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

template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

}; // namespace IDTO

#endif // UTILS_H
