#include "PZsparse.h"

namespace RAPTOR {
namespace Armour {

using namespace boost::multiprecision;

/*
Helper functions
*/

float getCenter(const Interval& a) {
    return (a.lower() + a.upper()) * 0.5;
}

float getRadius(const Interval& a) {
    return (a.upper() - a.lower()) * 0.5;
}

bool Monomial_sorter_degree(Monomial const& lhs, Monomial const& rhs) {
    return lhs.degree < rhs.degree;
}

/*
Initialization
*/

PZsparse::PZsparse() {
    center = 0;
    independent = 0;
}

PZsparse::PZsparse(const PZsparse& pz_inp) {
    center = pz_inp.center;
    polynomial = pz_inp.polynomial;
    independent = pz_inp.independent;
}

PZsparse::PZsparse(float center_inp) {
    center = center_inp;
    independent = 0;
}

PZsparse::PZsparse(float center_inp, float uncertainty_percent) {
    center = center_inp;
    independent = fabs(uncertainty_percent * center_inp);
}


PZsparse::PZsparse(float center_inp, Interval independent_inp) {
    center = center_inp + getCenter(independent_inp);
    independent = getRadius(independent_inp);
}

PZsparse::PZsparse(float center_inp, float* coeff_inp, uint32_t degree_inp[][NUM_FACTORS * 6], size_t num_monomials) {
    center = center_inp;

    polynomial.reserve(num_monomials);

    for (size_t i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = 0;

    simplify();
}

PZsparse::PZsparse(float center_inp, float* coeff_inp, uint32_t degree_inp[][NUM_FACTORS * 6], size_t num_monomials, Interval independent_inp) {
    center = center_inp + getCenter(independent_inp);

    polynomial.reserve(num_monomials);

    for (size_t i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    // assume independent_inp is centered at 0
    independent = getRadius(independent_inp);

    simplify();
}

/*
Internal functions
*/
void PZsparse::simplify() {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: simplify(): independent generator matrix has negative entry!");
    }

    sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);

    float reduce_amount = 0;

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    size_t i = 0;
    while (i < polynomial.size()) {
        size_t j = 0;
        cpp_int degree = polynomial[i].degree;

        for (j = i + 1; j < polynomial.size(); j++) {
            if (polynomial[j].degree != degree) {
                break;
            }

            polynomial[i].coeff += polynomial[j].coeff;
        }

        float temp = fabs(polynomial[i].coeff);
        if (temp <= SIMPLIFY_THRESHOLD) {
            reduce_amount += temp;
        }
        else {
            polynomial_new.emplace_back(polynomial[i]);
        }

        i = j;
    }

    if (!polynomial_new.empty()) {
        if (polynomial_new[0].degree == 0) {
            center += polynomial_new[0].coeff;
            polynomial.clear();
            polynomial.insert(polynomial.begin(), polynomial_new.begin() + 1, polynomial_new.end());
        }
        else {
            polynomial = polynomial_new;
        }
    }
    else {
        polynomial.clear();
    }

    if (reduce_amount != 0) {
        independent = independent + reduce_amount;
    }
}

// void PZsparse::reduce() {
//     if (independent < 0) {
//         throw std::runtime_error("PZsparse error: reduce(): independent generator matrix has negative entry!");
//     }

//     std::vector<Monomial> polynomial_new;
//     polynomial_new.reserve(polynomial.size());

//     for (auto it : polynomial) {
//         if (it.degree < max_hash_dependent_k_only) { // only dependent on k
//             polynomial_new.emplace_back(it.coeff, it.degree);
//         }
//         else {
//             independent += fabs(it.coeff);
//         }
//     }

//     polynomial = polynomial_new;
// }

Interval PZsparse::slice(const float factor[]) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: slice(): independent generator matrix has negative entry!");
    }
    float res_center = center;
    float res_radius = independent;

    uint32_t degreeArray[NUM_FACTORS * 6];

    for (auto it : polynomial) {
        float resTemp = it.coeff;

        if (it.degree < pow(MOVE_INC, NUM_FACTORS)) { // only dependent on k
            convertHashToDegree(degreeArray, it.degree);

            for (size_t j = 0; j < NUM_FACTORS; j++) {
                resTemp *= pow(factor[j], degreeArray[j]);
            }

            res_center += resTemp;
        }
        else { // this line should never be triggered if you run reduce first
            res_radius += fabs(resTemp);
        }
    }

    return Interval(res_center - res_radius, res_center + res_radius);
}

Interval PZsparse::slice(const Eigen::VectorXf& factor) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: slice(): independent generator matrix has negative entry!");
    }

    if (factor.size() != NUM_FACTORS) {
        throw std::runtime_error("PZsparse error: slice(): factor size does not match NUM_FACTORS!");
    }

    float res_center = center;
    float res_radius = independent;

    uint32_t degreeArray[NUM_FACTORS * 6];

    for (auto it : polynomial) {
        float resTemp = it.coeff;

        if (it.degree < pow(MOVE_INC, NUM_FACTORS)) { // only dependent on k
            convertHashToDegree(degreeArray, it.degree);

            for (size_t j = 0; j < NUM_FACTORS; j++) {
                resTemp *= pow(factor(j), degreeArray[j]);
            }

            res_center += resTemp;
        }
        else { // this line should never be triggered if you run reduce first
            res_radius += fabs(resTemp);
        }
    }

    return Interval(res_center - res_radius, res_center + res_radius);
}

void PZsparse::slice(float gradient[], const float factor[]) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: slice(): independent generator matrix has negative entry!");
    }

    std::memset(gradient, 0, NUM_FACTORS * sizeof(float));

    Eigen::Array<float, NUM_FACTORS, 1> resTemp;

    uint32_t degreeArray[NUM_FACTORS * 6];

    for (auto it : polynomial) {
        if (it.degree <= pow(MOVE_INC, NUM_FACTORS)) { // only dependent on k
            for (size_t k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff;
            }

            convertHashToDegree(degreeArray, it.degree);

            for (size_t j = 0; j < NUM_FACTORS; j++) {
                for (size_t k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k] = 0;
                        }
                        else {
                            resTemp[k] *= degreeArray[j] * pow(factor[j], degreeArray[j] - 1);
                        }
                    }
                    else {
                        resTemp[k] *= pow(factor[j], degreeArray[j]);
                    }
                }
            }

            for (size_t k = 0; k < NUM_FACTORS; k++) {
                gradient[k] += resTemp[k];
            }
        }
    }
}

Interval PZsparse::toInterval() {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: toInterval(): independent generator matrix has negative entry!");
    }

    float res_radius = independent;

    for (auto it : polynomial) {
        res_radius += fabs(it.coeff);
    }

    return Interval(center - res_radius, center + res_radius);
}

void convertHashToDegree(uint32_t degreeArray[], cpp_int degree) {
    for (size_t i = 0; i < NUM_FACTORS * 6; i++) {
        cpp_int res = degree % MOVE_INC;
        degreeArray[i] = res.convert_to<uint32_t>();    
        degree /= MOVE_INC;
    }

    return;
}

cpp_int convertDegreeToHash(const uint32_t degreeArray[]) {
    cpp_int degree = 0;
    cpp_int move_bit = 1;

    for (size_t i = 0; i < NUM_FACTORS * 6; i++) {
        if (degreeArray[i] > MOVE_INC) {
            throw std::runtime_error("degree can not be larger than " + std::to_string(MOVE_INC.convert_to<uint32_t>()) + "!");
        }
        degree += degreeArray[i] * move_bit;
        move_bit *= MOVE_INC;
    }

    return degree;
}

std::ostream& operator<<(std::ostream& os, const PZsparse& a) {
    os << a.center << " +...\n";

    uint32_t degreeArray[NUM_FACTORS * 6];

    for (auto it : a.polynomial) {
        os << '(' << it.coeff << ')';
        
        convertHashToDegree(degreeArray, it.degree);
        
        os << " * k^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j];
        }
        os << ") ";

        os << " * qde^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 1];
        }
        os << ") ";

        os << " * qdae^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 2];
        }
        os << ") ";

        os << " * qddae^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 3];
        }
        os << ") ";

        os << " * cosqe^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 4];
        }
        os << ") ";

        os << " * sinqe^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 5];
        }
        os << ") ";

        os << " +...\n";
    }

    os << "[ " << -a.independent << ", " << a.independent << " ]\n\n";

    return os;
}

std::ostream& operator<<(std::ostream& os, const Interval& a) {
    os << "[ " << a.lower() << ", " << a.upper() << "]\n";
    return os;
}

/*
Arithmetic
*/

PZsparse PZsparse::operator=(const float a) {
    center = a;
    polynomial.clear();
    independent = 0;
    return *this;
}

PZsparse PZsparse::operator=(const Interval& a) {
    center = getCenter(a);
    independent = getRadius(a);
    return *this;
}

PZsparse PZsparse::operator=(const PZsparse& a) {
    center = a.center;
    polynomial = a.polynomial;
    independent = a.independent;
    return *this;
}

PZsparse PZsparse::operator-() const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZsparse res;
    
    res.center = -center;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = independent;

    return res;
}

PZsparse PZsparse::operator+(const PZsparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator+(): independent generator matrix has negative entry!");
    }
    
    PZsparse res;

    res.center = center + a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    float reduce_amount = 0;
    size_t i = 0, j = 0;
 
    while (i < polynomial.size() && j < a.polynomial.size())
    {
        Monomial temp;
        if (polynomial[i].degree < a.polynomial[j].degree) {
            temp = polynomial[i];
            i++;
        }
        else if (polynomial[i].degree > a.polynomial[j].degree) {
            temp = a.polynomial[j];
            j++;
        }
        else {
            temp = Monomial(polynomial[i].coeff + a.polynomial[j].coeff, polynomial[i].degree);
            i++;
            j++;
        }

        if (fabs(temp.coeff) > SIMPLIFY_THRESHOLD) {
            res.polynomial.push_back(temp);
        }
        else {
            reduce_amount += fabs(temp.coeff);
        }
    }
 
    while (i < polynomial.size()) {
        if (fabs(polynomial[i].coeff) > SIMPLIFY_THRESHOLD) {
            res.polynomial.push_back(polynomial[i]);
        }
        else {
            reduce_amount += fabs(polynomial[i].coeff);
        }
        i++;
    }
 
    while (j < a.polynomial.size()) {
        if (fabs(a.polynomial[j].coeff) > SIMPLIFY_THRESHOLD) {
            res.polynomial.push_back(a.polynomial[j]);
        }
        else {
            reduce_amount += fabs(a.polynomial[j].coeff);
        }
        j++;
    }

    res.independent = independent + a.independent + reduce_amount;
    
    return res;
}

PZsparse PZsparse::operator+(const float a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator+(): independent generator matrix has negative entry!");
    }

    PZsparse res = *this;

    res.center += a;

    return res;
}

PZsparse operator+(const float a, const PZsparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZsparse error: operator+(): independent generator matrix has negative entry!");
    }

    PZsparse res = b;

    res.center += a;

    return res;
}

PZsparse PZsparse::operator+=(const PZsparse& a) {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator+=(): independent generator matrix has negative entry!");
    }
    
    center += a.center;

    polynomial.reserve(polynomial.size() + a.polynomial.size());

    // for (auto it : a.polynomial) {
    //     polynomial.push_back(it);
    // }

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size() + a.polynomial.size());

    float reduce_amount = 0;
    size_t i = 0, j = 0;
 
    while (i < polynomial.size() && j < a.polynomial.size())
    {
        Monomial temp;
        if (polynomial[i].degree < a.polynomial[j].degree) {
            temp = polynomial[i];
            i++;
        }
        else if (polynomial[i].degree > a.polynomial[j].degree) {
            temp = a.polynomial[j];
            j++;
        }
        else {
            temp = Monomial(polynomial[i].coeff + a.polynomial[j].coeff, polynomial[i].degree);
            i++;
            j++;
        }

        if (fabs(temp.coeff) > SIMPLIFY_THRESHOLD) {
            polynomial_new.push_back(temp);
        }
        else {
            reduce_amount += fabs(temp.coeff);
        }
    }
 
    while (i < polynomial.size()) {
        if (fabs(polynomial[i].coeff) > SIMPLIFY_THRESHOLD) {
            polynomial_new.push_back(polynomial[i]);
        }
        else {
            reduce_amount += fabs(polynomial[i].coeff);
        }
        i++;
    }
 
    while (j < a.polynomial.size()) {
        if (fabs(a.polynomial[j].coeff) > SIMPLIFY_THRESHOLD) {
            polynomial_new.push_back(a.polynomial[j]);
        }
        else {
            reduce_amount += fabs(a.polynomial[j].coeff);
        }
        j++;
    }

    polynomial = polynomial_new;

    independent += a.independent + reduce_amount;

    // simplify();
    
    return *this;
}

PZsparse PZsparse::operator-(const PZsparse& a) const{
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator-(): independent generator matrix has negative entry!");
    }    
    
    PZsparse res;

    res.center = center - a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    float reduce_amount = 0;
    size_t i = 0, j = 0;
 
    while (i < polynomial.size() && j < a.polynomial.size())
    {
        Monomial temp;
        if (polynomial[i].degree < a.polynomial[j].degree) {
            temp = polynomial[i];
            i++;
        }
        else if (polynomial[i].degree > a.polynomial[j].degree) {
            temp = a.polynomial[j];
            temp.coeff = -temp.coeff;
            j++;
        }
        else {
            temp = Monomial(polynomial[i].coeff - a.polynomial[j].coeff, polynomial[i].degree);
            i++;
            j++;
        }

        if (fabs(temp.coeff) > SIMPLIFY_THRESHOLD) {
            res.polynomial.push_back(temp);
        }
        else {
            reduce_amount += fabs(temp.coeff);
        }
    }
 
    while (i < polynomial.size()) {
        if (fabs(polynomial[i].coeff) > SIMPLIFY_THRESHOLD) {
            res.polynomial.push_back(polynomial[i]);
        }
        else {
            reduce_amount += fabs(polynomial[i].coeff);
        }
        i++;
    }
 
    while (j < a.polynomial.size()) {
        if (fabs(a.polynomial[j].coeff) > SIMPLIFY_THRESHOLD) {
            Monomial temp = a.polynomial[j];
            temp.coeff = -temp.coeff;
            res.polynomial.push_back(temp);
        }
        else {
            reduce_amount += fabs(a.polynomial[j].coeff);
        }
        j++;
    }

    res.independent = independent + a.independent + reduce_amount;
    
    return res;
}

PZsparse PZsparse::operator-(const float a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZsparse res = *this;

    res.center -= a;

    return res;
}

PZsparse operator-(const float a, const PZsparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZsparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZsparse res;

    res.center = a - b.center;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = b.independent;

    return res;
}

PZsparse PZsparse::operator*(const PZsparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator*(): independent generator matrix has negative entry!");
    }
    if (a.independent < 0) {
        throw std::runtime_error("PZsparse error: operator*(): independent generator matrix has negative entry!");
    }
    
    PZsparse res;

    // center * center
    res.center = center * a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size() + polynomial.size() * a.polynomial.size());

    // a.center * polynomial
    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff * a.center, it.degree);
    }

    // center * a.polynomial
    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(center * it.coeff, it.degree);
    }

    // polynomial * a.polynomial
    for (auto it1 : polynomial) {
        for (auto it2 : a.polynomial) {
            float multiply_coeff = it1.coeff * it2.coeff;

            multiply_coeff = it1.coeff * it2.coeff;

            // Do not have to check carry
            // if we already know the maximum degree in the polynomial
            res.polynomial.emplace_back(multiply_coeff, it1.degree + it2.degree);
        }
    }

    if ((!a.polynomial.empty()) && (!polynomial.empty())) {
        res.simplify();
    }

    // a.independent * (center + polynomial)
    float reduce_amount_2 = fabs(center);

    for (auto it : polynomial) {
        reduce_amount_2 += fabs(it.coeff);
    }

    reduce_amount_2 *= a.independent;
    
    // independent * (a.center + a.polynomial)
    float reduce_amount_3 = fabs(a.center);

    for (auto it : a.polynomial) {
        reduce_amount_3 += fabs(it.coeff);
    }
    
    reduce_amount_3 = independent * reduce_amount_3;

    // independent * a.independent + add reduced intervals
    float reduce_amount = reduce_amount_2 + reduce_amount_3;

    res.independent = independent * a.independent + reduce_amount;
    
    return res;
}

PZsparse PZsparse::operator*(const float a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZsparse res;

    res.center = center * a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = independent * fabs(a);

    return res;
}

PZsparse operator*(const float a, const PZsparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZsparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZsparse res;

    res.center = b.center * a;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = b.independent * fabs(a);

    return res;
}

PZsparse PZsparse::operator/(const float a) const {
    if (independent < 0) {
        throw std::runtime_error("PZsparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZsparse res;

    res.center = center / a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff / a, it.degree);
    }

    res.independent = independent / fabs(a);

    return res;
}

}; // namespace Armour
}; // namespace RAPTOR