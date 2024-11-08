#ifndef ARMOUR_HEADER_H
#define ARMOUR_HEADER_H

#include <cstdio>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <fstream>
#include <cstring>
#include <cassert>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <limits>

#include <Eigen/Dense>

#include <boost/numeric/interval.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace RAPTOR {
namespace Kinova {
namespace Armour {

// Boost intervals
namespace bn = boost::numeric;
namespace bi = bn::interval_lib;

using Interval = bn::interval<
    double, 
    bi::policies<
        bi::save_state<bi::rounded_transc_std<double>>,
        bi::checking_base<double>
    > 
>;

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // ARMOUR_HEADER_H