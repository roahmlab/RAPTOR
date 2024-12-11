#include <nanobind/nanobind.h>

#include "EndEffectorPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(end_effector_sysid_nanobind , m) {
    m.doc() = "nanobind end_effector_sysid plugin";

    nb::class_<EndEffectorPybindWrapper>(m, "EndEffectorPybindWrapper")
        .def(nb::init<const std::string,
                const std::string,
                const std::string, 
                const int,
                const bool>())
        .def("set_ipopt_parameters", &EndEffectorPybindWrapper::set_ipopt_parameters)
        .def("optimize", &EndEffectorPybindWrapper::optimize)
        .def("analyze_solution", &EndEffectorPybindWrapper::analyze_solution);
}
