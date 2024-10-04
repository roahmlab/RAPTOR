#include <nanobind/nanobind.h>

#include "DigitFillUnactuatedJointsPybindWrapper.h"

namespace nb = nanobind;

using namespace RAPTOR;
using namespace Digit;

NB_MODULE(DigitFillUnactuated_nanobind, m) {
    m.doc() = "nanobind DigitFillUnactuated_nanobind plugin";

    nb::class_<DigitFillUnactuatedJointsPybindWrapper>(m, "DigitFillUnactuatedJointsPybindWrapper")
        .def(nb::init<const char>())
        .def("fill_q", &DigitFillUnactuatedJointsPybindWrapper::fill_q)
        .def("fill_qv", &DigitFillUnactuatedJointsPybindWrapper::fill_qv)
        .def("fill_qva", &DigitFillUnactuatedJointsPybindWrapper::fill_qva);
}