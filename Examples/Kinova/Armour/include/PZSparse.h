#ifndef PZ_SPARSE_H
#define PZ_SPARSE_H

#include "ArmourHeaders.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

using namespace boost::multiprecision;

// make sure the number of monomials in a PZSparse is always smaller than this number
#define PZ_ORDER 20

// For now we only support 7 variables and the number is hardcoded unfortunately
// The degree of the robot has to be 7.
// If not, we need to change the following macros to support a different number of variables
#define NUM_FACTORS 7

// A specialized implementation of PZSparse
// so that more variables can be tracked
// k
// qe
// qde
// qdae
// qddae
// 5 * 7 = 35 variables in total
#define NUM_VARIABLES (NUM_FACTORS * 5)

// The order of the Taylor expansion for sin and cos
constexpr size_t SIN_TAYLOR_ORDER = 3;
constexpr size_t COS_TAYLOR_ORDER = 3;

const cpp_int MOVE_INC = 3 * std::max(SIN_TAYLOR_ORDER, COS_TAYLOR_ORDER) + 1;

double getCenter(const Interval& a);

double getRadius(const Interval& a);

class Monomial {
public:
    double coeff = 0;
    cpp_int degree = 0; // one hash unsigned integer for degrees for all factors, _K4_ _K3_ _K2_ _K1_

    Monomial() = default;
    Monomial(double coeff_inp, cpp_int degree_inp) : coeff(coeff_inp), degree(degree_inp) {};
};

class PZSparse {
public:
    // center
    double center = 0;

    // a vector of monomials for polynomial part
    std::vector<Monomial> polynomial;

    // intervals centered at 0 as independent generators
    // here we only store the upper bound of the intervals 
    // the actual interval is (-independent,independent)
    double independent = 0;

    /*
    Initialization
    */
    PZSparse();

    PZSparse(const PZSparse& pz_inp);

    PZSparse(double center_inp);

    PZSparse(double center_inp, double uncertainty_percent);

    PZSparse(const Interval a);

    PZSparse(double center_inp, Interval independent_inp);

    PZSparse(double center_inp, double* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials);

    PZSparse(double center_inp, double* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials, Interval independent_inp);

    ~PZSparse() = default;

    /*
    Internal functions
    */
    void simplify();

    void reduce();

    Interval slice(const double factor[]) const;

    Interval slice(const Eigen::VectorXd& factor) const;

    void slice(double gradient[], const double factor[]) const; // 1st-order gradient of slice of a 1-dim PZ

    void slice(Eigen::VectorXd& gradient, const double factor[]) const; // 1st-order gradient of slice of a 1-dim PZ

    Interval toInterval() const;

    /*
    Arithmetic
    */
    PZSparse operator=(const double a);
    
    PZSparse operator=(const Interval& a);

    PZSparse operator=(const PZSparse& a);

    PZSparse operator-() const;

    PZSparse operator-=(const PZSparse& a);

    PZSparse operator+(const PZSparse& a) const;

    PZSparse operator+(const double a) const;

    PZSparse operator+=(const PZSparse& a);

    PZSparse operator+=(const Interval& a);

    PZSparse operator-(const PZSparse& a) const;

    PZSparse operator-(const double a) const;

    PZSparse operator*(const PZSparse& a) const;

    PZSparse operator*=(const PZSparse& a);

    PZSparse operator*(const double a) const;

    PZSparse operator/(const double a) const;
    
    PZSparse operator/(const PZSparse& a) const;

    PZSparse operator/=(const PZSparse& a);

    bool operator<(const PZSparse& a) const;

    bool operator<(const double a) const;

    bool operator>(const PZSparse& a) const;

    bool operator>(const double a) const;

    bool operator==(const PZSparse& a) const;
};

typedef Eigen::Matrix<PZSparse, Eigen::Dynamic, Eigen::Dynamic> PZSparseMatrix;

void convertHashToDegree(uint32_t degreeArray[], const cpp_int degree);

cpp_int convertDegreeToHash(const uint32_t degreeArray[]);

std::ostream& operator<<(std::ostream& os, const PZSparse& a);

std::ostream& operator<<(std::ostream& os, const Interval& a);

/*
Arithmetic
*/

PZSparse operator+(const double a, const PZSparse& b);

PZSparse operator-(const double a, const PZSparse& b);

PZSparse operator*(const double a, const PZSparse& b);

PZSparse abs(const PZSparse& a);

PZSparse sqrt(const PZSparse& a);

bool isfinite(const PZSparse& a);

PZSparse sin(const PZSparse& a);

PZSparse cos(const PZSparse& a);

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

// define numeric_limits for PZSparse
namespace std {
    template<>
    class numeric_limits<RAPTOR::Kinova::Armour::PZSparse> {
    public:
        static constexpr bool is_specialized = true;
        static constexpr bool is_signed = true;
        static constexpr bool is_integer = false;
        static RAPTOR::Kinova::Armour::PZSparse min() noexcept {
            return RAPTOR::Kinova::Armour::PZSparse(numeric_limits<double>::min());
        }
        static RAPTOR::Kinova::Armour::PZSparse max() noexcept {
            return RAPTOR::Kinova::Armour::PZSparse(numeric_limits<double>::max());
        }
        static RAPTOR::Kinova::Armour::PZSparse lowest() noexcept {
            return RAPTOR::Kinova::Armour::PZSparse(numeric_limits<double>::lowest());
        }
        static RAPTOR::Kinova::Armour::PZSparse epsilon() noexcept {
            return RAPTOR::Kinova::Armour::PZSparse(numeric_limits<double>::epsilon());
        }
    };
}; // namespace std

#endif