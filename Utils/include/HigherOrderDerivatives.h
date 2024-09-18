#ifndef HIGHER_ORDER_DERIVATIVES_H
#define HIGHER_ORDER_DERIVATIVES_H

#include <cmath>
#include <stdexcept>

namespace RAPTOR {
namespace HigherOrderDerivatives {

inline float safeasin(const float x,
                       const bool throw_exception = false) {
    if (x > 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return M_PI / 2.0;
    } 
    else if (x < -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -M_PI / 2.0;
    } 
    return asinf(x);
}

inline float safedasindx(const float x,
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

inline float safeddasinddx(const float x,
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

inline float safedddasindddx(const float x,
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
    const float xSquare = x * x;
    return (1 + 2 * xSquare) / std::pow(1.0 - xSquare, 2.5);
}

inline float safeacos(const float x,
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
    return acosf(x);
}

inline float safedacosdx(const float x,
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

inline float safeddacosddx(const float x,
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

inline float safedddacosdddx(const float x,
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
    const float xSquare = x * x;
    return -(1 + 2 * xSquare) / std::pow(1.0 - xSquare, 2.5);
}

inline float safexSinx(const float x,
                        const float nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        float xSquare = x * x;
        float xFourth = xSquare * xSquare;
        return 1.0 + xSquare / 6.0 + 7.0 * xFourth / 360.0; // + O(x^6)
    }
    return x / sinf(x);
}

inline float safedxSinxdx(const float x,
                           float nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        float xSquare = x * x;
        float xThird = x * xSquare;
        float XFifth = xThird * xSquare;
        return x / 3.0 + 7.0 * xThird / 90.0 + 31.0 * XFifth / 2520.0; // + O(x^7)
    }
    const float sinx = sinf(x);
    return (sinx - x * cosf(x)) / (sinx * sinx);
}

inline float safeddxSinxddx(const float x,
                             const float nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        float xSquare = x * x;
        float xFourth = xSquare * xSquare;
        float xSixth = xFourth * xSquare;
        return 1.0 / 3.0 + 7.0 * xSquare / 30.0 + 31.0 * xFourth / 504.0; // + O(x^6)
    }
    const float sinx = sinf(x);
    const float cosx = cosf(x);
    const float sinxSquare = sinx * sinx;
    return (x * sinxSquare - 2 * cosx * sinx + 2 * x * cosx * cosx) / (sinxSquare * sinx);
}

inline float safedddxSinxdddx(const float x,
                               const float nearZeroThreshold = 1e-4) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        float xSquare = x * x;
        float xThird = x * xSquare;
        float XFifth = xThird * xSquare;
        return 7.0 * x / 15.0 + 31.0 * xThird / 126.0 + 127.0 * XFifth / 1800.0; // + O(x^7)
    }
    const float sinx = sinf(x);
    const float cosx = cosf(x);
    const float sinxSquare = sinx * sinx;
    return x * cosx / sinxSquare +
           6 * cosx * cosx / (sinx * sinxSquare) -
           6 * x * cosx / (sinxSquare * sinxSquare) + 
           3 / sinx;
}

}; // namespace HigherOrderDerivatives
}; // namespace RAPTOR

#endif // HIGHER_ORDER_DERIVATIVES_H