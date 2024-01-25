#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "KinovaPybindWrapper.h"

namespace py = pybind11;

using namespace IDTO;
using namespace Kinova;

PYBIND11_MODULE(oracle_pybind, m) {
    m.doc() = "pybind11 oracle_pybind plugin";

    py::class_<KinovaPybindWrapper>(m, "KinovaPybindWrapper")
        .def(py::init<const std::string>())
        .def("set_obstacles", &KinovaPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &KinovaPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &KinovaPybindWrapper::set_trajectory_parameters)
        .def("set_buffer", &KinovaPybindWrapper::set_buffer)
        .def("set_target", &KinovaPybindWrapper::set_target)
        .def("optimize", &KinovaPybindWrapper::optimize);
}