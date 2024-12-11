#include <nanobind/nanobind.h>

#include "EndEffectorIdentificationPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(end_effector_sysid_nanobind , m) {
    m.doc() = "nanobind end effector identification plugin";

    nb::class_<EndEffectorIdentificationPybindWrapper>(m, "EndEffectorIdentificationPybindWrapper")
        .def(nb::init<const std::string, const std::string, const std::string, const int, const bool>())
        .def("set_ipopt_parameters", &EndEffectorIdentificationPybindWrapper::set_ipopt_parameters)
        .def("optimize", &EndEffectorIdentificationPybindWrapper::optimize);
}
