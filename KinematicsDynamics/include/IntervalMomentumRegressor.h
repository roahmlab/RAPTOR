#ifndef INTERVAL_MOMENTUM_REGRESSOR_H
#define INTERVAL_MOMENTUM_REGRESSOR_H

#include "IntervalRegressorInverseDynamics.h"

namespace RAPTOR {

// Compute inverse dynamics using tau = Y * phi,
// where Y is the n x 10*n regressor matrix and 
// phi is the vector of 10*n dynamic parameters (inertia, com, mass for each of the link).
// phi is constant and directly loaded from the robot. 
// The gradient of Y will also be computed.
class IntervalMomentumRegressor : public IntervalRegressorInverseDynamics {
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
    IntervalMomentumRegressor() = default;

    // Constructor
    IntervalMomentumRegressor(const Model& model_input, 
                              const std::shared_ptr<TrajectoryData>& trajPtr_input,
                              const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
                              Eigen::VectorXi jtype_input = Eigen::VectorXi(0)) :
        IntervalRegressorInverseDynamics(model_input, trajPtr_input, sensor_noise_input, jtype_input) {
        Y_CTv.resize(N * modelPtr_->nv, numParams);
        pY_CTv_pz.resize(trajPtr_->varLength);
        for (int i = 0; i < trajPtr_->varLength; i++) {
            pY_CTv_pz(i).resize(N * modelPtr_->nv, numParams);
        }
    };

    // Destructor
    ~IntervalMomentumRegressor() = default;

    // class methods:
    virtual void compute(const VecXd& z,
                         bool compute_derivatives = true) final override;

    // class members:

    // this has been defined in IntervalRegressorInverseDynamics 
    // and stores the regressor for system momentum H(q) * v
    // MatXInt Y;
    // Eigen::Array<MatXInt, 1, Eigen::Dynamic> pY_pz;

    // this stores the regressor for C^T(q, v) * v
    // which is needed to compute the time derivative of the system momentum
    MatXInt Y_CTv;
    Eigen::Array<MatXInt, 1, Eigen::Dynamic> pY_CTv_pz;
};

}; // namespace RAPTOR

#endif // INTERVAL_MOMENTUM_REGRESSOR_H