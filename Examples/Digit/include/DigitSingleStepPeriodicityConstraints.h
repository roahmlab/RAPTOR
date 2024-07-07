
#ifndef DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include "Constraints.h"
#include "DigitConstrainedInverseDynamics.h"
#include "SurfaceContactConstraints.h"

#include <memory>
#include <string>

namespace IDTO {
namespace Digit {    

// ... ==> ... 
const std::string JOINT_MAP[NUM_INDEPENDENT_JOINTS][2] = {
    {"left_hip_roll",   "right_hip_roll"}, 
    {"left_hip_yaw",    "right_hip_yaw"},
    {"left_hip_pitch",  "right_hip_pitch"}, 
    {"left_knee",       "right_knee"},
    {"left_toe_A",      "right_toe_A"}, 
    {"left_toe_B",      "right_toe_B"},
    {"right_hip_roll",  "left_hip_roll"}, 
    {"right_hip_yaw",   "left_hip_yaw"},
    {"right_hip_pitch", "left_hip_pitch"}, 
    {"right_knee",      "left_knee"},
    {"right_toe_A",     "left_toe_A"}, 
    {"right_toe_B",     "left_toe_B"},
};

class DigitSingleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitSingleStepPeriodicityConstraints() = default;

    // Constructor
    DigitSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                          std::shared_ptr<ConstrainedInverseDynamics> dcidPtr_input,
                                          const frictionParams& fp_input);

    // Destructor
    ~DigitSingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

    void compute_bounds() override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<ConstrainedInverseDynamics> dcidPtr_;

    frictionParams fp;

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

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
