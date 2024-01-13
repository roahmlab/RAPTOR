#ifndef PYBIND_BASE_HPP
#define PYBIND_BASE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace IDTO {

namespace py = pybind11;

class PybindBase {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    PybindBase();

    // Destructor
    ~PybindBase();

    // Class methods
    virtual py::tuple optimize() = 0;
};

}; // namespace IDTO

#endif // PYBIND_BASE_HPP