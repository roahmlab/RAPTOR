#include <nanobind/nanobind.h>

#include "SafePayloadExcitingPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

NB_MODULE(safe_payload_exciting_nanobind, m) {
    m.doc() = "nanobind armour_nanobind plugin";

    nb::class_<SafePayloadExcitingPybindWrapper>(m, "SafePayloadExcitingPybindWrapper")
        .def(nb::init<const std::string, const std::string, const bool>())
        .def("set_obstacles", &SafePayloadExcitingPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &SafePayloadExcitingPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &SafePayloadExcitingPybindWrapper::set_trajectory_parameters)
        .def("set_endeffector_inertial_parameters", &SafePayloadExcitingPybindWrapper::set_endeffector_inertial_parameters)
        .def("optimize", &SafePayloadExcitingPybindWrapper::optimize)
        .def("analyze_solution", &SafePayloadExcitingPybindWrapper::analyze_solution);
}