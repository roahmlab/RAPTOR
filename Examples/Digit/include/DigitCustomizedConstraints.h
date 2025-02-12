#ifndef DIGIT_CUSTOMIZED_CONSTRAINTS_H
#define DIGIT_CUSTOMIZED_CONSTRAINTS_H

#include "Constraints.h"
#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

namespace RAPTOR {
namespace Digit {

typedef struct GaitParameters_ {
    double eps_torso_angle = Utils::deg2rad(1.5); // 1.5 degrees
    double swingfoot_midstep_z_des = 0.15; // meters
    double swingfoot_begin_x_des = -0.22; // meters (negative if left stance, positive if right stance)
    double swingfoot_begin_y_des = 0.00; // meters
    double swingfoot_end_x_des = -0.22; // meters (negative if left stance, positive if right stance)
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
                               std::shared_ptr<Trajectories>& trajPtr_input,
                               std::shared_ptr<DigitDynamicsConstraints>& dcPtr_input,
                               const GaitParameters& gp_input);

    // Destructor
    ~DigitCustomizedConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override;

        // print violation info
    virtual void print_violation_info() override;

    // class variables:
    GaitParameters gp;

    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<DigitDynamicsConstraints> ddcPtr_;

    std::unique_ptr<Model> modelPtr_;

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_;

        // swing foot information
    Model::JointIndex swingfoot_id = 0;

    Transform swingfoot_endT;

    Transform leftfoot_endT;
    Transform rightfoot_endT;

        // updated in compute()
    MatX q;
    MatX swingfoot_xyzrpy;

    Eigen::Array<MatX, 1, Eigen::Dynamic> pq_pz;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pswingfoot_xyzrpy_pz;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;

        // constraints
    VecX g1, g1_lb, g1_ub;
    VecX g2, g2_lb, g2_ub;
    VecX g3, g3_lb, g3_ub;
    VecX g4, g4_lb, g4_ub;
    VecX g5, g5_lb, g5_ub;
    VecX g6, g6_lb, g6_ub;
    VecX g7, g7_lb, g7_ub;
    VecX g8, g8_lb, g8_ub;
    VecX g9, g9_lb, g9_ub;
    VecX g10, g10_lb, g10_ub;
    VecX g11, g11_lb, g11_ub;
};

} // namespace Digit
} // namespace RAPTOR

#endif // DIGIT_CUSTOMIZED_CONSTRAINTS_H
