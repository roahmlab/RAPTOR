
#ifndef G1_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define G1_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include "Constraints.h"
#include "G1ConstrainedInverseDynamics.h"
#include "RectangleSurfaceContactConstraints.h"

#include <memory>
#include <string>

namespace RAPTOR {
namespace G1 {    

// ... ==> ... 
const std::string JOINT_MAP[NUM_INDEPENDENT_JOINTS][2] = {
    {"left_hip_pitch_joint",    "right_hip_pitch_joint"},
    {"left_hip_roll_joint",     "right_hip_roll_joint"},
    {"left_hip_yaw_joint",      "right_hip_yaw_joint"},
    {"left_knee_joint",         "right_knee_joint"},
    {"left_ankle_pitch_joint",  "right_ankle_pitch_joint"},
    {"left_ankle_roll_joint",   "right_ankle_roll_joint"},
    {"right_hip_pitch_joint",   "left_hip_pitch_joint"},
    {"right_hip_roll_joint",    "left_hip_roll_joint"},
    {"right_hip_yaw_joint",     "left_hip_yaw_joint"},
    {"right_knee_joint",        "left_knee_joint"},
    {"right_ankle_pitch_joint", "left_ankle_pitch_joint"},
    {"right_ankle_roll_joint",  "left_ankle_roll_joint"}
};

class G1SingleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    G1SingleStepPeriodicityConstraints() = default;

    // Constructor
    G1SingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                       std::shared_ptr<G1ConstrainedInverseDynamics> dcidPtr_input,
                                       const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~G1SingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<G1ConstrainedInverseDynamics> dcidPtr_;

    rectangleContactSurfaceParams fp;

        // intermediate variables updated in compute()
    int joint_id1[NUM_INDEPENDENT_JOINTS];
    int joint_id2[NUM_INDEPENDENT_JOINTS];

    MatX prnea_pq;
    MatX prnea_pv;
    MatX prnea_pa;
    VecX zeroVec;

    MatX pv_plus_pz;
    MatX plambda_pz;

    VecX g1;
    VecX g2;
    VecX g3;
    VecX g4;
    VecX g5;

    MatX pg1_pz;
    MatX pg2_pz;
    MatX pg3_pz;
    MatX pg4_pz;
    MatX pg5_pz;
};

}; // namespace G1
}; // namespace RAPTOR

#endif // G1_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
