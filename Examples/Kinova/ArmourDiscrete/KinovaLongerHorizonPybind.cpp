#include <nanobind/nanobind.h>

#include "KinovaLongerHorizonPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(kinova_longer_nanobind, m) {
    m.doc() = "nanobind kinova_longer_nanobind plugin";

    nb::class_<KinovaLongerHorizonPybindWrapper>(m, "KinovaLongerHorizonPybindWrapper")
        .def(nb::init<const std::string, const bool>())
        .def("set_obstacles", &KinovaLongerHorizonPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &KinovaLongerHorizonPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &KinovaLongerHorizonPybindWrapper::set_trajectory_parameters)
        .def("set_buffer", &KinovaLongerHorizonPybindWrapper::set_buffer)
        .def("optimize", &KinovaLongerHorizonPybindWrapper::optimize)
        .def("analyze_solution", &KinovaLongerHorizonPybindWrapper::analyze_solution);
}