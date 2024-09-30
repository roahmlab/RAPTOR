#ifndef PZ_SPARSE_H
#define PZ_SPARSE_H

#include "ArmourHeaders.h"

namespace RAPTOR {
namespace Armour {

using namespace boost::multiprecision;

// monomials with a coefficient smaller than this number will be reduced
#define SIMPLIFY_THRESHOLD 1e-4

// For now we only support 7 variables and the number is hardcoded unfortunately
// The degree of the robot has to be 7.
// If not, we need to change the following macros to support a different number of variables
#define NUM_FACTORS 7

// A specialized implementation of PZsparse
// so that more variables can be tracked
// k
// qe
// cosqe
// sinqe
// qde
// qdae
// qddae
// 7 * 7 = 49 variables in total
#define NUM_VARIABLES (NUM_FACTORS * 7)

const cpp_int MOVE_INC = 7;

// const cpp_int max_hash_dependent_k_only = ((cpp_int)1 << (cpp_int)(2 * NUM_FACTORS));

// const cpp_int max_hash_dependent_k_links_only = ((cpp_int)1 << (cpp_int)(5 * NUM_FACTORS));
// const cpp_int dependent_k_mask = max_hash_dependent_k_only - 1;

float getCenter(const Interval& a);

float getRadius(const Interval& a);

float getCenter(const Interval& a);

float getRadius(const Interval& a);

class Monomial {
public:
    float coeff = 0;
    cpp_int degree = 0; // one hash unsigned integer for degrees for all factors, _K4_ _K3_ _K2_ _K1_

    Monomial() = default;
    Monomial(float coeff_inp, cpp_int degree_inp) : coeff(coeff_inp), degree(degree_inp) {};
};

class PZsparse {
public:
    // center
    float center = 0;

    // a vector of monomials for polynomial part
    std::vector<Monomial> polynomial;

    // intervals centered at 0 as independent generators
    // here we only store the upper bound of the intervals 
    // the actual interval is (-independent,independent)
    float independent = 0;

    /*
    Initialization
    */
    PZsparse();

    PZsparse(const PZsparse& pz_inp);

    PZsparse(float center_inp);

    PZsparse(float center_inp, float uncertainty_percent);

    PZsparse(float center_inp, Interval independent_inp);

    PZsparse(float center_inp, float* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials);

    PZsparse(float center_inp, float* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials, Interval independent_inp);

    ~PZsparse() = default;

    /*
    Internal functions
    */
    void simplify();

    void reduce();

    Interval slice(const float factor[]) const;

    Interval slice(const Eigen::VectorXf& factor) const;

    void slice(float gradient[], const float factor[]) const; // 1st-order gradient of slice of a 1-dim PZ

    void slice(Eigen::VectorXf& gradient, const float factor[]) const; // 1st-order gradient of slice of a 1-dim PZ

    Interval toInterval();

    /*
    Arithmetic
    */
    PZsparse operator=(const float a);
    
    PZsparse operator=(const Interval& a);

    PZsparse operator=(const PZsparse& a);

    PZsparse operator-() const;

    PZsparse operator+(const PZsparse& a) const;

    PZsparse operator+(const float a) const;

    PZsparse operator+=(const PZsparse& a);

    PZsparse operator-(const PZsparse& a) const;

    PZsparse operator-(const float a) const;

    PZsparse operator*(const PZsparse& a) const;

    PZsparse operator*(const float a) const;

    PZsparse operator/(const float a) const;
};

typedef Eigen::Array<PZsparse, Eigen::Dynamic, Eigen::Dynamic> PZsparseArray;

void convertHashToDegree(uint32_t degreeArray[], const cpp_int degree);

cpp_int convertDegreeToHash(const uint32_t degreeArray[]);

std::ostream& operator<<(std::ostream& os, const PZsparse& a);

std::ostream& operator<<(std::ostream& os, const Interval& a);

/*
Arithmetic
*/

PZsparse operator+(const float a, const PZsparse& b);

PZsparse operator-(const float a, const PZsparse& b);

PZsparse operator*(const float a, const PZsparse& b);

}; // namespace Armour
}; // namespace RAPTOR

#endif