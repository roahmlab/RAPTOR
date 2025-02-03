#ifndef INTERVAL_MOMENTUM_REGRESSOR_H
#define INTERVAL_MOMENTUM_REGRESSOR_H

#include "MomentumRegressor.h"
#include "IntervalRegressorInverseDynamics.h"

namespace RAPTOR {

// Compute inverse dynamics using tau = Y * phi,
// where Y is the n x 10*n regressor matrix and 
// phi is the vector of 10*n dynamic parameters (inertia, com, mass for each of the link).
// phi is constant and directly loaded from the robot. 
// The gradient of Y will also be computed.
class IntervalMomentumRegressor : public MomentumRegressor {
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
                              Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~IntervalMomentumRegressor() = default;

    // class methods:
    virtual void compute(const VecXd& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override;

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

    MatXInt Y_CTv;
    Eigen::Array<MatXInt, 1, Eigen::Dynamic> pY_CTv_pz;

    Eigen::Array<VecXInt, 1, Eigen::Dynamic> tau;
    Eigen::Array<MatXInt, 1, Eigen::Dynamic> ptau_pz;

    int NB = 0;

    std::shared_ptr<MomentumRegressor> mrPtr_ = nullptr;
};

}; // namespace RAPTOR

#endif // INTERVAL_MOMENTUM_REGRESSOR_H