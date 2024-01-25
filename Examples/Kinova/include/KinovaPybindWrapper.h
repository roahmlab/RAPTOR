#ifndef KINOVA_PYBIND_WRAPPER_H
#define KINOVA_PYBIND_WRAPPER_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "KinovaOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace IDTO {
namespace Kinova {

namespace py = pybind11;

class KinovaPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    KinovaPybindWrapper() = default;

    KinovaPybindWrapper(const std::string urdf_filename);

    // Destructor
    ~KinovaPybindWrapper() = default;

    // Class methods
    void set_obstacles(const int num_obstacles_inp,
                       const py::array_t<double> obstacles_inp);

    void set_ipopt_parameters(const double tol,
                              const double obj_scaling_factor,
                              const double max_wall_time, 
                              const int print_level,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    void set_trajectory_parameters(const py::array_t<double> q0_inp,
                                   const py::array_t<double> qd0_inp,
                                   const py::array_t<double> qdd0_inp,
                                   const double duration_inp);

    void set_buffer(const py::array_t<double> joint_limits_buffer_inp,
                    const py::array_t<double> velocity_limits_buffer_inp,
                    const py::array_t<double> torque_limits_buffer_inp);

    void set_target(const py::array_t<double> q_des_inp,
                    const double tplan_inp);

    py::tuple optimize();

    // Class members
    // robot model
    Model model;
    Eigen::VectorXi jtype;

    // obstacle information
    int num_obstacles = 0;
    Eigen::Array<Vec3, 1, Eigen::Dynamic> zonotopeCenters;
    Eigen::Array<MatX, 1, Eigen::Dynamic> zonotopeGenerators;

    // trajectory information
    ArmourTrajectoryParameters atp;
    double T = 1;
    int N = 16;
    int degree = 5;
    VecX qdes;
    double tplan = 0;
    int tplan_n = 0;

    // buffer information
    VecX joint_limits_buffer;
    VecX velocity_limits_buffer;
    VecX torque_limits_buffer;
    
    SmartPtr<KinovaOptimizer> mynlp;
    SmartPtr<IpoptApplication> app;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_ipopt_parameters_check = false;
    bool set_trajectory_parameters_check = false;
    bool set_buffer_check = false;
    bool set_target_check = false;
};

}; // namespace Kinova
}; // namespace IDTO

#endif // KINOVA_PYBIND_WRAPPER_H
