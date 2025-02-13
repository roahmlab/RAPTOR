#include <nanobind/nanobind.h>

#include "EndEffectorIdentificationMomentumPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Kinova;

NB_MODULE(end_effector_sysid_momentum_nanobind , m) {
    m.doc() = "nanobind end effector identification plugin";

    nb::class_<EndEffectorIdentificationMomentumPybindWrapper>(m, "EndEffectorIdentificationMomentumPybindWrapper")
        .def(nb::init<const std::string, 
                      const nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>, 
                      const int, 
                      const std::string, 
                      const bool>())
        .def("set_ipopt_parameters", &EndEffectorIdentificationMomentumPybindWrapper::set_ipopt_parameters)
        .def("add_trajectory_file", &EndEffectorIdentificationMomentumPybindWrapper::add_trajectory_file)
        .def("optimize", &EndEffectorIdentificationMomentumPybindWrapper::optimize)
        .def("reset", &EndEffectorIdentificationMomentumPybindWrapper::reset);
}
