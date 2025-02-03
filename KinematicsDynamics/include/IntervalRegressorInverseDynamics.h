#ifndef INTERVAL_REGRESSOR_INVERSEDYNAMICS_H
#define INTERVAL_REGRESSOR_INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/regressor.hpp"

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp> 

#include "RegressorInverseDynamics.h"
#include "CustomizedInverseDynamics.h"
#include "Spatial.h"
#include "TrajectoryData.h"

#include <cmath>
#include <iostream> 
#include <memory>
#include <cstdio>
#include <cstdlib>

namespace RAPTOR {

// Boost interval arithmetic
namespace bn = boost::numeric;
namespace bi = bn::interval_lib;

using Interval = bn::interval<
    double, 
    bi::policies<
        bi::save_state<bi::rounded_transc_std<double>>,
        bi::checking_base<double>
    > 
>;

namespace IntervalHelper {
double getCenter(const Interval& x);
double getRadius(const Interval& x);
Interval makeErrorInterval(const double error, 
                           const SensorNoiseInfo::SensorNoiseType type, 
                           const double value);
}; // namespace IntervalHelper

// Compute inverse dynamics using tau = Y * phi,
// where Y is the n x 10*n regressor matrix and 
// phi is the vector of 10*n dynamic parameters (inertia, com, mass for each of the link).
// phi is constant and directly loaded from the robot. 
// The gradient of Y will also be computed.
class IntervalRegressorInverseDynamics : public RegressorInverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecXInt = Eigen::Vector<Interval, Eigen::Dynamic>;
    using MatXInt = Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec3Int = Eigen::Vector<Interval, 3>;
    using Mat3Int = Eigen::Matrix<Interval, 3, 3>;
    using Vec6Int = Eigen::Vector<Interval, 6>;
    using Mat6Int = Eigen::Matrix<Interval, 6, 6>;
    using MatRegressorInt = Eigen::Matrix<Interval, 6, 10>;
    using VecXd = Eigen::VectorXd;
    using MatXd = Eigen::MatrixXd;
    using Vec6d = Eigen::Vector<double, 6>;
    using Mat6d = Eigen::Matrix<double, 6, 6>;

    // Constructor
    IntervalRegressorInverseDynamics() = default;

    // Constructor
    IntervalRegressorInverseDynamics(const Model& model_input, 
                                     const std::shared_ptr<TrajectoryData>& trajPtr_input,
                                     Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~IntervalRegressorInverseDynamics() = default;

    // class methods:
    virtual void compute(const VecXd& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;
                                           
    // class members:
    std::shared_ptr<Model> modelPtr_ = nullptr;
    std::shared_ptr<Data> dataPtr_ = nullptr;

    std::shared_ptr<TrajectoryData> trajPtr_ = nullptr;

    int N = 0; // number of time instances in tspan

    Eigen::VectorXi jtype;
    Eigen::Array<Mat6d, 1, Eigen::Dynamic> Xtree;
    Vec6d a_grav;
    
    VecXd phi;

    int numParams = 0;

    MatXInt Y;
    Eigen::Array<MatXInt, 1, Eigen::Dynamic> pY_pz;

    Eigen::Array<VecXInt, 1, Eigen::Dynamic> tau;
    Eigen::Array<MatXInt, 1, Eigen::Dynamic> ptau_pz;

    int NB = 0;

    std::shared_ptr<RegressorInverseDynamics> ridPtr_ = nullptr;
};

// muliply operator between Eigen interval matrices and vectors
Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> intervalMatrixMultiply(
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& B);

Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> intervalDoubleMatrixMultiply(
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& B);

Eigen::Matrix<Interval, 6, 6> intervalDouble66MatrixMultiply(
    const Eigen::Matrix<Interval, 6, 6>& A,
    const Eigen::Matrix<double, 6, 6>& B);

Eigen::Matrix<Interval, 6, 1> intervalDouble61VectorMultiply(
    const Eigen::Matrix<Interval, 6, 6>& A,
    const Eigen::Matrix<double, 6, 1>& B);

}; // namespace RAPTOR

#endif // INTERVAL_REGRESSOR_INVERSEDYNAMICS_H