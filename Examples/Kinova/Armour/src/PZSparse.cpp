#include "PZSparse.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

using namespace boost::multiprecision;

/*
Helper functions
*/

double getCenter(const Interval& a) {
    return (a.lower() + a.upper()) * 0.5;
}

double getRadius(const Interval& a) {
    return (a.upper() - a.lower()) * 0.5;
}

bool Monomial_sorter_degree(Monomial const& lhs, Monomial const& rhs) {
    return lhs.degree < rhs.degree;
}

bool Monomial_sorter_coeff(Monomial const& lhs, Monomial const& rhs) {
    return fabs(lhs.coeff) > fabs(rhs.coeff);
}

/*
Initialization
*/

PZSparse::PZSparse() {
    center = 0;
    independent = 0;
}

PZSparse::PZSparse(const PZSparse& pz_inp) {
    center = pz_inp.center;
    polynomial = pz_inp.polynomial;
    independent = pz_inp.independent;
}

PZSparse::PZSparse(double center_inp) {
    center = center_inp;
    independent = 0;
}

PZSparse::PZSparse(double center_inp, double uncertainty_percent) {
    center = center_inp;
    independent = fabs(uncertainty_percent * center_inp);
}

PZSparse::PZSparse(const Interval a) {
    center = getCenter(a);
    independent = getRadius(a);
}

PZSparse::PZSparse(double center_inp, Interval independent_inp) {
    center = center_inp + getCenter(independent_inp);
    independent = getRadius(independent_inp);
}

PZSparse::PZSparse(double center_inp, double* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials) {
    center = center_inp;

    polynomial.reserve(num_monomials);

    for (size_t i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = 0;

    simplify();
}

PZSparse::PZSparse(double center_inp, double* coeff_inp, uint32_t degree_inp[][NUM_VARIABLES], size_t num_monomials, Interval independent_inp) {
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
void PZSparse::simplify() {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: simplify(): independent generator matrix has negative entry!");
    }

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    // sort by degree to combine monomials with the same degree
    std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);

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

        polynomial_new.emplace_back(polynomial[i]);

        i = j;
    }

    // assign the new polynomial
    if (!polynomial_new.empty()) {
        // if the first monomial is a constant, add it to the center
        if (polynomial_new[0].degree == 0) {
            center += polynomial_new[0].coeff;
            polynomial.clear();
            polynomial.insert(polynomial.begin(), polynomial_new.begin() + 1, polynomial_new.end());
        }
        else {
            polynomial = polynomial_new;
        }
    }
    else { // all monomials are eliminated
        polynomial.clear();
    }

    // reduce the number of monomials within PZ_ORDER
    if (polynomial.size() > PZ_ORDER) {
        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_coeff);

        double reduce_amount = 0;
        for (size_t i = PZ_ORDER; i < polynomial.size(); i++) {
            reduce_amount += fabs(polynomial[i].coeff);
        }

        polynomial.erase(polynomial.begin() + PZ_ORDER, polynomial.end());

        if (reduce_amount != 0) {
            independent += reduce_amount;
        }

        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);
    }
}

void PZSparse::reduce() {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: reduce(): independent generator matrix has negative entry!");
    }

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    for (auto it : polynomial) {
        if (it.degree < pow(MOVE_INC, NUM_FACTORS)) { // only dependent on k
            polynomial_new.emplace_back(it.coeff, it.degree);
        }
        else {
            independent += fabs(it.coeff);
        }
    }

    polynomial = polynomial_new;
}

Interval PZSparse::slice(const double factor[]) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: slice(): independent generator matrix has negative entry!");
    }
    double res_center = center;
    double res_radius = independent;

    uint32_t degreeArray[NUM_VARIABLES];

    for (auto it : polynomial) {
        double resTemp = it.coeff;

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

Interval PZSparse::slice(const Eigen::VectorXd& factor) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: slice(): independent generator matrix has negative entry!");
    }

    if (factor.size() != NUM_FACTORS) {
        throw std::runtime_error("PZSparse error: slice(): factor size does not match NUM_FACTORS!");
    }

    double res_center = center;
    double res_radius = independent;

    uint32_t degreeArray[NUM_VARIABLES];

    for (auto it : polynomial) {
        double resTemp = it.coeff;

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

void PZSparse::slice(double gradient[], const double factor[]) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: slice(): independent generator matrix has negative entry!");
    }

    std::memset(gradient, 0, NUM_FACTORS * sizeof(double));

    Eigen::VectorXd resTemp(NUM_FACTORS);

    uint32_t degreeArray[NUM_VARIABLES];

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

void PZSparse::slice(Eigen::VectorXd& gradient, const double factor[]) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: slice(): independent generator matrix has negative entry!");
    }

    if (gradient.size() != NUM_FACTORS) {
        throw std::invalid_argument("PZSparse error: slice(): gradient size does not match NUM_FACTORS!");
    }

    gradient.setZero();

    Eigen::VectorXd resTemp(NUM_FACTORS);

    uint32_t degreeArray[NUM_VARIABLES];

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

Interval PZSparse::toInterval() const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: toInterval(): independent generator matrix has negative entry!");
    }

    double res_radius = independent;

    for (auto it : polynomial) {
        res_radius += fabs(it.coeff);
    }

    return Interval(center - res_radius, center + res_radius);
}

void convertHashToDegree(uint32_t degreeArray[], cpp_int degree) {
    for (size_t i = 0; i < NUM_VARIABLES; i++) {
        cpp_int res = degree % MOVE_INC;
        degreeArray[i] = res.convert_to<uint32_t>();    
        degree /= MOVE_INC;
    }

    return;
}

cpp_int convertDegreeToHash(const uint32_t degreeArray[]) {
    cpp_int degree = 0;
    cpp_int move_bit = 1;

    for (size_t i = 0; i < NUM_VARIABLES; i++) {
        if (degreeArray[i] > MOVE_INC) {
            throw std::runtime_error("degree can not be larger than " + std::to_string(MOVE_INC.convert_to<uint32_t>()) + "!");
        }
        degree += degreeArray[i] * move_bit;
        move_bit *= MOVE_INC;
    }

    return degree;
}

std::ostream& operator<<(std::ostream& os, const PZSparse& a) {
    os << a.center << " +...\n";

    uint32_t degreeArray[NUM_VARIABLES];

    for (auto it : a.polynomial) {
        os << '(' << it.coeff << ')';
        
        convertHashToDegree(degreeArray, it.degree);
        
        os << " * k^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j];
        }
        os << ") ";

        os << " * qe^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 1];
        }
        os << ") ";

        os << " * qde^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 2];
        }
        os << ") ";

        os << " * qdae^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 3];
        }
        os << ") ";

        os << " * qddae^(";
        for (size_t j = 0; j < NUM_FACTORS; j++) {
            os << degreeArray[j + NUM_FACTORS * 4];
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

PZSparse PZSparse::operator=(const double a) {
    center = a;
    polynomial.clear();
    independent = 0;
    return *this;
}

PZSparse PZSparse::operator=(const Interval& a) {
    center = getCenter(a);
    independent = getRadius(a);
    return *this;
}

PZSparse PZSparse::operator=(const PZSparse& a) {
    center = a.center;
    polynomial = a.polynomial;
    independent = a.independent;
    return *this;
}

PZSparse PZSparse::operator-() const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZSparse res;
    
    res.center = -center;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = independent;

    return res;
}

PZSparse PZSparse::operator-=(const PZSparse& a) {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator-=(): independent generator matrix has negative entry!");
    }
    
    center -= a.center;

    polynomial.reserve(polynomial.size() + a.polynomial.size());

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size() + a.polynomial.size());

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

        if (temp.coeff != 0) {
            polynomial_new.push_back(temp);
        }
    }
 
    while (i < polynomial.size()) {
        polynomial_new.push_back(polynomial[i]);
        i++;
    }
 
    while (j < a.polynomial.size()) {
        Monomial temp = a.polynomial[j];
        temp.coeff = -temp.coeff;
        polynomial_new.push_back(temp);
        j++;
    }

    polynomial = polynomial_new;

    // reduce the number of monomials within PZ_ORDER
    if (polynomial.size() > PZ_ORDER) {
        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_coeff);

        double reduce_amount = 0;
        for (size_t i = PZ_ORDER; i < polynomial.size(); i++) {
            reduce_amount += fabs(polynomial[i].coeff);
        }

        polynomial.erase(polynomial.begin() + PZ_ORDER, polynomial.end());

        if (reduce_amount != 0) {
            independent += reduce_amount;
        }

        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);
    }
    
    return *this;
}

PZSparse PZSparse::operator+(const PZSparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator+(): independent generator matrix has negative entry!");
    }
    
    PZSparse res;

    res.center = center + a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

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

        if (temp.coeff != 0) {
            res.polynomial.push_back(temp);
        }
    }
 
    while (i < polynomial.size()) {
        res.polynomial.push_back(polynomial[i]);
        i++;
    }
 
    while (j < a.polynomial.size()) {
        res.polynomial.push_back(a.polynomial[j]);
        j++;
    }

    res.independent = independent + a.independent;

    // reduce the number of monomials within PZ_ORDER
    if (res.polynomial.size() > PZ_ORDER) {
        std::sort(res.polynomial.begin(), res.polynomial.end(), Monomial_sorter_coeff);

        double reduce_amount = 0;
        for (size_t i = PZ_ORDER; i < res.polynomial.size(); i++) {
            reduce_amount += fabs(res.polynomial[i].coeff);
        }

        res.polynomial.erase(res.polynomial.begin() + PZ_ORDER, res.polynomial.end());

        if (reduce_amount != 0) {
            res.independent += reduce_amount;
        }

        std::sort(res.polynomial.begin(), res.polynomial.end(), Monomial_sorter_degree);
    }
    
    return res;
}

PZSparse PZSparse::operator+(const double a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator+(): independent generator matrix has negative entry!");
    }

    PZSparse res = *this;

    res.center += a;

    return res;
}

PZSparse operator+(const double a, const PZSparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZSparse error: operator+(): independent generator matrix has negative entry!");
    }

    PZSparse res = b;

    res.center += a;

    return res;
}

PZSparse PZSparse::operator+=(const PZSparse& a) {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator+=(): independent generator matrix has negative entry!");
    }
    
    center += a.center;

    polynomial.reserve(polynomial.size() + a.polynomial.size());

    std::vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size() + a.polynomial.size());

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

        if (temp.coeff != 0) {
            polynomial_new.push_back(temp);
        }
    }
 
    while (i < polynomial.size()) {
        polynomial_new.push_back(polynomial[i]);
        i++;
    }
 
    while (j < a.polynomial.size()) {
        polynomial_new.push_back(a.polynomial[j]);
        j++;
    }

    polynomial = polynomial_new;
    independent += a.independent;
    
    // reduce the number of monomials within PZ_ORDER
    if (polynomial.size() > PZ_ORDER) {
        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_coeff);

        double reduce_amount = 0;
        for (size_t i = PZ_ORDER; i < polynomial.size(); i++) {
            reduce_amount += fabs(polynomial[i].coeff);
        }

        polynomial.erase(polynomial.begin() + PZ_ORDER, polynomial.end());

        if (reduce_amount != 0) {
            independent += reduce_amount;
        }

        std::sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);
    }
    
    return *this;
}

PZSparse PZSparse::operator+=(const Interval& a) {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator+=(): independent generator matrix has negative entry!");
    }

    center += getCenter(a);
    independent += getRadius(a);

    return *this;
}

PZSparse PZSparse::operator-(const PZSparse& a) const{
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator-(): independent generator matrix has negative entry!");
    }    
    
    PZSparse res;

    res.center = center - a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

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

        if (temp.coeff != 0) {
            res.polynomial.push_back(temp);
        }
    }
 
    while (i < polynomial.size()) {
        res.polynomial.push_back(polynomial[i]);
        i++;
    }
 
    while (j < a.polynomial.size()) {
        Monomial temp = a.polynomial[j];
        temp.coeff = -temp.coeff;
        res.polynomial.push_back(temp);
        j++;
    }

    res.independent = independent + a.independent;

    // reduce the number of monomials within PZ_ORDER
    if (res.polynomial.size() > PZ_ORDER) {
        std::sort(res.polynomial.begin(), res.polynomial.end(), Monomial_sorter_coeff);

        double reduce_amount = 0;
        for (size_t i = PZ_ORDER; i < res.polynomial.size(); i++) {
            reduce_amount += fabs(res.polynomial[i].coeff);
        }

        res.polynomial.erase(res.polynomial.begin() + PZ_ORDER, res.polynomial.end());

        if (reduce_amount != 0) {
            res.independent += reduce_amount;
        }

        std::sort(res.polynomial.begin(), res.polynomial.end(), Monomial_sorter_degree);
    }
    
    return res;
}

PZSparse PZSparse::operator-(const double a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZSparse res = *this;

    res.center -= a;

    return res;
}

PZSparse operator-(const double a, const PZSparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZSparse error: operator-(): independent generator matrix has negative entry!");
    }

    PZSparse res;

    res.center = a - b.center;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = b.independent;

    return res;
}

PZSparse PZSparse::operator*(const PZSparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator*(): independent generator matrix has negative entry!");
    }
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: operator*(): independent generator matrix has negative entry!");
    }
    
    PZSparse res;

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
            double multiply_coeff = it1.coeff * it2.coeff;

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
    double reduce_amount_2 = fabs(center);

    for (auto it : polynomial) {
        reduce_amount_2 += fabs(it.coeff);
    }

    reduce_amount_2 *= a.independent;
    
    // independent * (a.center + a.polynomial)
    double reduce_amount_3 = fabs(a.center);

    for (auto it : a.polynomial) {
        reduce_amount_3 += fabs(it.coeff);
    }
    
    reduce_amount_3 = independent * reduce_amount_3;

    // independent * a.independent + add reduced intervals
    double reduce_amount = reduce_amount_2 + reduce_amount_3;

    res.independent = independent * a.independent + reduce_amount;
    
    return res;
}

PZSparse PZSparse::operator*=(const PZSparse& a) {
    PZSparse res = *this * a;
    *this = res;
    return *this;
}

PZSparse PZSparse::operator*(const double a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZSparse res;

    res.center = center * a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = independent * fabs(a);

    return res;
}

PZSparse operator*(const double a, const PZSparse& b) {
    if (b.independent < 0) {
        throw std::runtime_error("PZSparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZSparse res;

    res.center = b.center * a;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = b.independent * fabs(a);

    return res;
}

PZSparse PZSparse::operator/(const double a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator*(): independent generator matrix has negative entry!");
    }

    PZSparse res;

    res.center = center / a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff / a, it.degree);
    }

    res.independent = independent / fabs(a);

    return res;
}

PZSparse PZSparse::operator/(const PZSparse& a) const {
    // a is just a scalar
    if (a.polynomial.size() == 0 && a.independent == 0) {
        return *this / a.center;
    }

    std::cerr << a << std::endl;
    throw std::runtime_error("PZSparse error: operator/(): not supported!");
    return PZSparse(0);
}

PZSparse PZSparse::operator/=(const PZSparse& a) {
    // a is just a scalar
    if (a.polynomial.size() == 0 && a.independent == 0) {
        *this = *this / a.center;
        return *this;
    }

    std::cerr << a << std::endl;
    throw std::runtime_error("PZSparse error: operator/=(): not supported!");
    return *this;
}

bool PZSparse::operator<(const PZSparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator<(): independent generator matrix has negative entry!");
    }
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: operator<(): independent generator matrix has negative entry!");
    }

    return center < a.center;
}

bool PZSparse::operator<(const double a) const {
    return center < a;
}

bool PZSparse::operator>(const PZSparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator>(): independent generator matrix has negative entry!");
    }
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: operator>(): independent generator matrix has negative entry!");
    }

    return center > a.center;
}

bool PZSparse::operator>(const double a) const {
    return center > a;
} 

bool PZSparse::operator==(const PZSparse& a) const {
    if (independent < 0) {
        throw std::runtime_error("PZSparse error: operator==(): independent generator matrix has negative entry!");
    }
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: operator==(): independent generator matrix has negative entry!");
    }
    
    if (center != a.center) {
        return false;
    }

    if (independent != a.independent) {
        return false;
    }

    if (polynomial.size() != a.polynomial.size()) {
        return false;
    }

    for (size_t i = 0; i < polynomial.size(); i++) {
        if (polynomial[i].coeff != a.polynomial[i].coeff) {
            return false;
        }
        if (polynomial[i].degree != a.polynomial[i].degree) {
            return false;
        }
    }

    return true;
}

PZSparse abs(const PZSparse& a) {
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: abs(): independent generator matrix has negative entry!");
    }

    PZSparse res;
    res.center = fabs(a.center);

    res.polynomial.reserve(a.polynomial.size());

    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(fabs(it.coeff), it.degree);
    }

    res.independent = a.independent;

    return res;
}

PZSparse sqrt(const PZSparse& a) {
    // a is just a scalar
    if (a.polynomial.size() == 0 && a.independent == 0) {
        return PZSparse(std::sqrt(a.center));
    }

    std::cerr << a << std::endl;
    throw std::runtime_error("PZSparse error: sqrt(): not supported!");
    return PZSparse(0);
}

bool isfinite(const PZSparse& a) {
    if (a.independent < 0) {
        throw std::runtime_error("PZSparse error: isfinite(): independent generator matrix has negative entry!");
    }

    if (!std::isfinite(a.center)) {
        return false;
    }

    for (auto it : a.polynomial) {
        if (!std::isfinite(it.coeff)) {
            return false;
        }
    }

    if (!std::isfinite(a.independent)) {
        return false;
    }

    return true;
}

PZSparse sin(const PZSparse& a) {
    if (a.polynomial.size() == 0 && a.independent == 0) {
        return PZSparse(std::sin(a.center));
    }
    
    PZSparse res(std::sin(a.center));

    // Taylor expansion
    const PZSparse multiplier = a - a.center;
    PZSparse power(1.0);
    for (size_t i = 0; i <= SIN_TAYLOR_ORDER; i++) {
        double coeff = 0;
        switch (i % 4) {
            case 0:
                coeff = std::cos(a.center); // 1st order derivative of sin
                break;
            case 1:
                coeff = -std::sin(a.center); // 2nd order derivative of sin
                break;
            case 2:
                coeff = -std::cos(a.center); // 3rd order derivative of sin
                break;
            case 3:
                coeff = std::sin(a.center); // 4th order derivative of sin
                break;
        }
        power *= multiplier / (double)(i + 1);

        if (i < SIN_TAYLOR_ORDER) {
            res += power * coeff;
        }
        else { // Lagrange remainder
            res += power.toInterval() * coeff;
        }
    }

    return res;
}

PZSparse cos(const PZSparse& a) {
    if (a.polynomial.size() == 0 && a.independent == 0) {
        return PZSparse(std::cos(a.center));
    }
    
    PZSparse res(std::cos(a.center));

    // Taylor expansion
    const PZSparse multiplier = a - a.center;
    PZSparse power(1.0);
    for (size_t i = 0; i <= SIN_TAYLOR_ORDER; i++) {
        double coeff = 0;
        switch (i % 4) {
            case 0:
                coeff = -std::sin(a.center); // 1st order derivative of cos
                break;
            case 1:
                coeff = -std::cos(a.center); // 2nd order derivative of cos
                break;
            case 2:
                coeff = std::sin(a.center); // 3rd order derivative of cos
                break;
            case 3:
                coeff = std::cos(a.center); // 4th order derivative of cos
                break;
        }
        power *= multiplier / (double)(i + 1);

        if (i < SIN_TAYLOR_ORDER) {
            res += power * coeff;
        }
        else { // Lagrange remainder
            res += power.toInterval() * coeff;
        }
    }

    return res;
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR