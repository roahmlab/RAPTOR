#include <nanobind/nanobind.h>

#include "EndEffectorIdentificationPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(end_effector_sysid_nanobind , m) {
    m.doc() = "nanobind end effector identification plugin";

    nb::class_<EndEffectorIdentificationPybindWrapper>(m, "EndEffectorIdentificationPybindWrapper")
        .def(nb::init<const std::string, 
                      const nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>,
                      const std::string, 
                      const bool>())
        .def("set_ipopt_parameters", &EndEffectorIdentificationPybindWrapper::set_ipopt_parameters)
        .def("add_trajectory_file", &EndEffectorIdentificationPybindWrapper::add_trajectory_file)
        .def("optimize", &EndEffectorIdentificationPybindWrapper::optimize)
        .def("reset", &EndEffectorIdentificationPybindWrapper::reset);
}
