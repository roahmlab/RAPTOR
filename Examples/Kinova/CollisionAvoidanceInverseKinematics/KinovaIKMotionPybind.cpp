#include <nanobind/nanobind.h>

#include "KinovaIKMotionPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(KinovaIKMotion_nanobind, m) {
    m.doc() = "nanobind KinovaIKMotion_nanobind plugin";

    nb::class_<KinovaIKMotionPybindWrapper>(m, "KinovaIKMotionPybindWrapper")
        .def(nb::init<const std::string, const bool>())
        .def("set_desired_endeffector_transforms", &KinovaIKMotionPybindWrapper::set_desired_endeffector_transforms)
        .def("set_obstacles", &KinovaIKMotionPybindWrapper::set_obstacles)
        .def("set_ipopt_parameters", &KinovaIKMotionPybindWrapper::set_ipopt_parameters)
        .def("solve", &KinovaIKMotionPybindWrapper::solve);
}