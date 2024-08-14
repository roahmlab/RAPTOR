
#ifndef DIGIT_MODIFIED_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define DIGIT_MODIFIED_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include "Constraints.h"
#include "DigitModifiedConstrainedInverseDynamics.h"
#include "RectangleSurfaceContactConstraints.h"

#include <memory>
#include <string>

namespace RAPTOR {
namespace DigitModified {    

// ... ==> ... 
const std::string JOINT_MAP[NUM_INDEPENDENT_JOINTS][2] = {
    {"left_hip_roll",   "right_hip_roll"}, 
    {"left_hip_yaw",    "right_hip_yaw"},
    {"left_hip_pitch",  "right_hip_pitch"}, 
    {"left_knee",       "right_knee"},
    {"left_tarsus",     "right_tarsus"},
    {"left_toe_pitch",  "right_toe_pitch"}, 
    {"left_toe_roll",   "right_toe_roll"},
    {"right_hip_roll",  "left_hip_roll"}, 
    {"right_hip_yaw",   "left_hip_yaw"},
    {"right_hip_pitch", "left_hip_pitch"}, 
    {"right_knee",      "left_knee"},
    {"right_tarsus",    "left_tarsus"},
    {"right_toe_pitch", "left_toe_pitch"}, 
    {"right_toe_roll",  "left_toe_roll"},
};

class DigitModifiedSingleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitModifiedSingleStepPeriodicityConstraints() = default;

    // Constructor
    DigitModifiedSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                  std::shared_ptr<DigitModifiedConstrainedInverseDynamics> dcidPtr_input,
                                                  const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~DigitModifiedSingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

    void compute_bounds() override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<DigitModifiedConstrainedInverseDynamics> dcidPtr_;

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

}; // namespace DigitModified
}; // namespace RAPTOR

#endif // DIGIT_MODIFIED_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
