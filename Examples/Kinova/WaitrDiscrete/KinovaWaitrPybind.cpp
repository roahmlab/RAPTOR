#include <nanobind/nanobind.h>

#include "KinovaWaitrPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(oracle_waitr_nanobind, m) {
    m.doc() = "nanobind oracle_waitr_nanobind plugin";

    nb::class_<KinovaWaitrPybindWrapper>(m, "KinovaWaitrPybindWrapper")
        .def(nb::init<const std::string>())
        .def("set_obstacles", &KinovaWaitrPybindWrapper::set_obstacles)
        .def("set_contact_surface_parameters", &KinovaWaitrPybindWrapper::set_contact_surface_parameters)
        .def("set_ipopt_parameters", &KinovaWaitrPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &KinovaWaitrPybindWrapper::set_trajectory_parameters)
        .def("set_buffer", &KinovaWaitrPybindWrapper::set_buffer)
        .def("set_target", &KinovaWaitrPybindWrapper::set_target)
        .def("optimize", &KinovaWaitrPybindWrapper::optimize)
        .def("analyze_solution", &KinovaWaitrPybindWrapper::analyze_solution);
}