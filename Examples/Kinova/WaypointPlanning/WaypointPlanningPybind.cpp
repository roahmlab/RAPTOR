#include <nanobind/nanobind.h>

#include "WaypointPlanningPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(KinovaHLP_nanobind, m) {
    m.doc() = "nanobind KinovaHLP_nanobind plugin";

    nb::class_<WaypointPlanningPybindWrapper>(m, "WaypointPlanningPybindWrapper")
        .def(nb::init<const std::string>())
        .def("set_obstacles", &WaypointPlanningPybindWrapper::set_obstacles)
        .def("set_start_goal", &WaypointPlanningPybindWrapper::set_start_goal)
        .def("plan", &WaypointPlanningPybindWrapper::plan);
}