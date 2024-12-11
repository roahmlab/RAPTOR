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

class EndEffectorPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_float = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_float = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
   

    // Constructor
    EndEffectorPybindWrapper() = default;

    EndEffectorPybindWrapper(const std::string urdf_filename,
                             const std::string trajectory_filename,
                             const std::string friction_parameters_filename,
                             const int H_input,
                             const bool display_info);

    // Destructor
    ~EndEffectorPybindWrapper() = default;

    // Class methods
    void set_ipopt_parameters(const double tol,
                          const double max_wall_time,
                          const int print_level,
                          const int max_iter,
                          const std::string mu_strategy,
                          const std::string linear_solver,
                          const bool gradient_check);

    Eigen::VectorXd x_to_theta(Eigen::VectorXd& z);
    Eigen::MatrixXd matrixSquareRoot(Eigen::MatrixXd& matrix);
    double BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B);
    Eigen::MatrixXd theta_to_LMI(Eigen::VectorXd& theta) ;

    nb::tuple optimize();

    double analyze_solution();

    // Class members
    // robot model
    Model model;

    // trajectory information
    const std::string trajectory_filename;
    const int H;
    const VecX offset;
     Eigen::VectorXd theta_sol;

    SmartPtr<EndEffectorParametersIdentification> mynlp;
    SmartPtr<IpoptApplication> app;

    // Flags to check if the parameters are set
    bool set_ipopt_parameters_check = false;
    bool has_optimized = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // ENDEFFECTOR_PYBIND_WRAPPER_H
