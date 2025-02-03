#include <nanobind/nanobind.h>

#include "ArmourPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

NB_MODULE(armour_nanobind, m) {
    m.doc() = "nanobind armour_nanobind plugin";

    nb::class_<ArmourPybindWrapper>(m, "ArmourPybindWrapper")
        .def(nb::init<const std::string, const std::string, const bool>())
        .def("set_obstacles", &ArmourPybindWrapper::set_obstacles)
        .def("set_endeffector_inertial_parameters", &ArmourPybindWrapper::set_endeffector_inertial_parameters)
        .def("set_ipopt_parameters", &ArmourPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &ArmourPybindWrapper::set_trajectory_parameters)
        .def("optimize", &ArmourPybindWrapper::optimize)
        .def("analyze_solution", &ArmourPybindWrapper::analyze_solution);
}