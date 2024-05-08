
#ifndef DIGIT_CUSTOMIZED_CONSTRAINTS_WITHOBSTACLES_H
#define DIGIT_CUSTOMIZED_CONSTRAINTS_WITHOBSTACLES_H

#include "Constraints.h"
#include "FourierCurves.h"
#include "ZonotopeCollisionAvoidance.h"
#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

namespace IDTO {
namespace Digit {

const double swingfoot_sphere_radius = 0.2350 * 0.5; // this is the length of Digit's foot
const double eps_torso_angle = 0.0524; // 3 degrees

class DigitCustomizedConstraintsWithObstacles : public Constraints {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitCustomizedConstraintsWithObstacles() = default;

    // Constructor
    DigitCustomizedConstraintsWithObstacles(const Model& model_input,
                                            const Eigen::VectorXi& jtype_input,
                                            std::shared_ptr<Trajectories>& trajPtr_input,
                                            std::shared_ptr<DynamicsConstraints>& dcPtr_input,
                                            const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                                            const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input,
                                            const VecX& q_act0_input,
                                            const VecX& q_act_d0_input);

    // Destructor
    ~DigitCustomizedConstraintsWithObstacles() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<Trajectories> trajPtr_ = nullptr;
    std::shared_ptr<DynamicsConstraints> dcPtr_ = nullptr;

    std::unique_ptr<Model> modelPtr_ = nullptr;

    std::unique_ptr<ForwardKinematicsHighOrderDerivative> fkhofPtr_ = nullptr;

    std::shared_ptr<ZonotopeCollisionAvoidance> collisionAvoidancePtr_ = nullptr;

        // jtype copy
    Eigen::VectorXi jtype;

        // the joint index of the joint we want to constrain
    Model::JointIndex swingfoot_id = 0;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform swingfoot_endT;

        // updated in compute()
    Transform jointT;
    MatX jointTJ;

    MatX q;
    MatX swingfoot_xyzrpy;

    Eigen::Array<MatX, 1, Eigen::Dynamic> pq_pz;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pswingfoot_xyzrpy_pz;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;

        // initial condition
    VecX q_act0;
    VecX q_act_d0;
};

} // namespace Digit
} // namespace IDTO

#endif // DIGIT_CUSTOMIZED_CONSTRAINTS_WITHOBSTACLES_H
