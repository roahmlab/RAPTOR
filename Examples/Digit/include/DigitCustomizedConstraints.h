
#ifndef DIGIT_CUSTOMIZED_CONSTRAINTS_H
#define DIGIT_CUSTOMIZED_CONSTRAINTS_H

#include "Constraints.h"
#include "FourierCurves.h"
#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

namespace IDTO {
namespace Digit {

typedef struct GaitParameters_ {
    double eps_torso_angle = 0.0524; // 3 degrees
    double swingfoot_midstep_z_des = 0.15; // meters
    double swingfoot_begin_x_des = -0.22; // meters
    double swingfoot_begin_y_des = 0.00; // meters
    double swingfoot_end_x_des = -0.22; // meters
    double swingfoot_end_y_des = 0.00; // meters
} GaitParameters;

class DigitCustomizedConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitCustomizedConstraints() = default;

    // Constructor
    DigitCustomizedConstraints(const Model& model_input,
                               const Eigen::VectorXi& jtype_input,
                               std::shared_ptr<Trajectories>& trajPtr_input,
                               std::shared_ptr<DynamicsConstraints>& dcPtr_input,
                               const GaitParameters& gp_input);

    // Destructor
    ~DigitCustomizedConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    GaitParameters gp;

    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<DynamicsConstraints> dcPtr_;

    std::unique_ptr<Model> modelPtr_;

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_;

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
};

} // namespace Digit
} // namespace IDTO

#endif // DIGIT_CUSTOMIZED_CONSTRAINTS_H
