#include <nanobind/nanobind.h>

#include "DualArmourPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

NB_MODULE(dualarmour_nanobind, m) {
    m.doc() = "nanobind dualarmour_nanobind plugin";

    nb::class_<DualArmourPybindWrapper>(m, "DualArmourPybindWrapper")
        .def(nb::init<const std::string, const std::string, const std::string, const std::string, const bool>())
        .def("set_obstacles", &DualArmourPybindWrapper::set_obstacles)
        .def("set_endeffector_inertial_parameters_robot1", &DualArmourPybindWrapper::set_endeffector_inertial_parameters_robot1)
        .def("set_endeffector_inertial_parameters_robot2", &DualArmourPybindWrapper::set_endeffector_inertial_parameters_robot2)
        .def("set_ipopt_parameters", &DualArmourPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &DualArmourPybindWrapper::set_trajectory_parameters)
        .def("optimize", &DualArmourPybindWrapper::optimize)
        .def("analyze_solution", &DualArmourPybindWrapper::analyze_solution);
}