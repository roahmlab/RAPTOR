#include <nanobind/nanobind.h>

#include "KinovaPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(kinova_nanobind, m) {
    m.doc() = "nanobind kinova_nanobind plugin";

    nb::class_<KinovaPybindWrapper>(m, "KinovaPybindWrapper")
        .def(nb::init<const std::string, const bool>())
        .def("set_obstacles", &KinovaPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &KinovaPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &KinovaPybindWrapper::set_trajectory_parameters)
        .def("set_buffer", &KinovaPybindWrapper::set_buffer)
        .def("set_target", &KinovaPybindWrapper::set_target)
        .def("optimize", &KinovaPybindWrapper::optimize)
        .def("analyze_solution", &KinovaPybindWrapper::analyze_solution);
}