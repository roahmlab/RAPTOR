#ifndef ENDEFFECTOR_PYBIND_WRAPPER_H
#define ENDEFFECTOR_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "EndEffectorParametersIdentification.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;

class EndEffectorIdentificationPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_float = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_float = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
   
    // Constructor
    EndEffectorIdentificationPybindWrapper() = default;

    EndEffectorIdentificationPybindWrapper(const std::string urdf_filename,
                                           const std::string trajectory_filename,
                                           const std::string friction_parameters_filename,
                                           const int H_input,
                                           const bool display_info);

    // Destructor
    ~EndEffectorIdentificationPybindWrapper() = default;

    // class methods
    void set_ipopt_parameters(const double tol,
                              const double max_wall_time,
                              const int print_level,
                              const int max_iter,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    nb::tuple optimize();

    // class members
        // robot model
    Model model;
    VecX offset;

        // trajectory information
    std::string trajectory_filename;

        // system identification settings
    int H = 10; // forward integration horizon
    int downsample_rate = 1; // 1 means no downsample

    SmartPtr<EndEffectorParametersIdentification> mynlp;
    SmartPtr<IpoptApplication> app;

    SensorNoiseInfo sensor_noise;

    // Flags to check if the parameters are set
    bool set_ipopt_parameters_check = false;
    bool has_optimized = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // ENDEFFECTOR_PYBIND_WRAPPER_H
