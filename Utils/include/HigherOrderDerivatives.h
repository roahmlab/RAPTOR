#ifndef HIGHER_ORDER_DERIVATIVES_H
#define HIGHER_ORDER_DERIVATIVES_H

#include <cmath>
#include <stdexcept>

namespace RAPTOR {
namespace HigherOrderDerivatives {

inline double safeasin(const double x,
                       const bool throw_exception = false) {
    if (x > 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return M_PI_2;
    } 
    else if (x < -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -M_PI_2;
    } 
    return asin(x);
}

inline double safedasindx(const double x,
                          const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return 1e10; // a very large positive number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    }    
    return 1.0 / std::sqrt(1.0 - x * x);             
}

inline double safeddasinddx(const double x,
                            const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return 1e10; // a very large positive number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    return x / std::pow(1.0 - x * x, 1.5);
}

inline double safedddasindddx(const double x,
                              const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return 1e10; // a very large positive number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    const double xSquare = x * x;
    return (1 + 2 * xSquare) / std::pow(1.0 - xSquare, 2.5);
}

inline double safeacos(const double x,
                       const bool throw_exception = false) {
    if (x > 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return 0.0;
    } 
    else if (x < -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return M_PI;
    } 
    return acos(x);
}

inline double safedacosdx(const double x,
                          const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return -1e10; // a very large negative number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return 1e10; // a very large positive number
    } 
    return -1.0 / std::sqrt(1.0 - x * x);
}

inline double safeddacosddx(const double x,
                            const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return 1e10; // a very large positive number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    return -x / std::pow(1.0 - x * x, 1.5);
}

inline double safedddacosdddx(const double x,
                              const bool throw_exception = false) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return -1e10; // a very large negative number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return 1e10; // a very large positive number
    } 
    const double xSquare = x * x;
    return -(1 + 2 * xSquare) / std::pow(1.0 - xSquare, 2.5);
}

inline double safexSinx(const double x,
                        const double nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xFourth = xSquare * xSquare;
        return 1.0 + xSquare / 6.0 + 7.0 * xFourth / 360.0; // + O(x^6)
    }
    return x / sin(x);
}

inline double safedxSinxdx(const double x,
                           double nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xThird = x * xSquare;
        double XFifth = xThird * xSquare;
        return x / 3.0 + 7.0 * xThird / 90.0 + 31.0 * XFifth / 2520.0; // + O(x^7)
    }
    const double sinx = sin(x);
    return (sinx - x * cos(x)) / (sinx * sinx);
}

inline double safeddxSinxddx(const double x,
                             const double nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xFourth = xSquare * xSquare;
        double xSixth = xFourth * xSquare;
        return 1.0 / 3.0 + 7.0 * xSquare / 30.0 + 31.0 * xFourth / 504.0; // + O(x^6)
    }
    const double sinx = sin(x);
    const double cosx = cos(x);
    const double sinxSquare = sinx * sinx;
    return (x * sinxSquare - 2 * cosx * sinx + 2 * x * cosx * cosx) / (sinxSquare * sinx);
}

inline double safedddxSinxdddx(const double x,
                               const double nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xThird = x * xSquare;
        double XFifth = xThird * xSquare;
        return 7.0 * x / 15.0 + 31.0 * xThird / 126.0 + 127.0 * XFifth / 1800.0; // + O(x^7)
    }
    const double sinx = sin(x);
    const double cosx = cos(x);
    const double sinxSquare = sinx * sinx;
    return x * cosx / sinxSquare +
           6 * cosx * cosx / (sinx * sinxSquare) -
           6 * x * cosx / (sinxSquare * sinxSquare) + 
           3 / sinx;
}

}; // namespace HigherOrderDerivatives
}; // namespace RAPTOR

#endif // HIGHER_ORDER_DERIVATIVES_H