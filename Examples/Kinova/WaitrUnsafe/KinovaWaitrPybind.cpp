#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "KinovaWaitrPybindWrapper.h"

namespace py = pybind11;

using namespace RAPTOR;
using namespace Kinova;

PYBIND11_MODULE(oracle_waitr_pybind, m) {
    m.doc() = "pybind11 oracle_waitr_pybind plugin";

    py::class_<KinovaWaitrPybindWrapper>(m, "KinovaWaitrPybindWrapper")
        .def(py::init<const std::string>())
        .def("set_obstacles", &KinovaWaitrPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &KinovaWaitrPybindWrapper::set_ipopt_parameters)
        .def("set_trajectory_parameters", &KinovaWaitrPybindWrapper::set_trajectory_parameters)
        .def("set_buffer", &KinovaWaitrPybindWrapper::set_buffer)
        .def("set_target", &KinovaWaitrPybindWrapper::set_target)
        .def("optimize", &KinovaWaitrPybindWrapper::optimize);
}