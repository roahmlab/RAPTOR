#ifndef ENDEFFECTOR_ID_MOMENTUM_PYBIND_WRAPPER_H
#define ENDEFFECTOR_ID_MOMENTUM_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "EndEffectorParametersIdentificationMomentum.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;

class EndEffectorIdentificationMomentumPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using Vec10 = Eigen::Vector<double, 10>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
   
    // Constructor
    EndEffectorIdentificationMomentumPybindWrapper() = default;

    EndEffectorIdentificationMomentumPybindWrapper(const std::string urdf_filename,
                                                   const nb_1d_double friction_parameters_input,
                                                   const int H_input,
                                                   const std::string time_format_string,
                                                   const bool display_info);

    // Destructor
    ~EndEffectorIdentificationMomentumPybindWrapper() = default;

    // class methods
    void set_ipopt_parameters(const double tol,
                              const double max_wall_time,
                              const int print_level,
                              const int max_iter,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    void add_trajectory_file(const std::string filename_input);

    nb::tuple optimize();

    void reset();

    // class members
        // robot model
    Model model;
    VecX offset;

        // trajectory information
    TimeFormat time_format;

        // system identification settings
    int H = 10; // forward integration horizon
    int downsample_rate = 1; // 1 means no downsample

    SmartPtr<EndEffectorParametersIdentificationMomentum> mynlp;
    SmartPtr<IpoptApplication> app;

    SensorNoiseInfo sensor_noise;

    // Flags to check if the parameters are set
    bool set_ipopt_parameters_check = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // ENDEFFECTOR_ID_MOMENTUM_PYBIND_WRAPPER_H
