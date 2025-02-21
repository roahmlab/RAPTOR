
#ifndef Biped_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define Biped_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include "Constraints.h"
#include "BipedConstrainedInverseDynamics.h"
#include "RectangleSurfaceContactConstraints.h"

#include <memory>
#include <string>

namespace RAPTOR {
namespace Biped {    

// ... ==> ... 
const std::string JOINT_MAP[NUM_INDEPENDENT_JOINTS][2] = {
    {"left_leg_1",  "right_leg_1"},
    {"left_leg_2",  "right_leg_2"},
    {"left_leg_3",  "right_leg_3"},
    {"right_leg_1", "left_leg_1"},
    {"right_leg_2", "left_leg_2"},
    {"right_leg_3", "left_leg_3"},
    {"shoulders",   "shoulders"},
    {"head",        "head"},
    {"left_arm_1",  "right_arm_1"},
    {"left_arm_2",  "right_arm_2"},
    {"left_arm_3",  "right_arm_3"},
    {"right_arm_1", "left_arm_1"},
    {"right_arm_2", "left_arm_2"},
    {"right_arm_3", "left_arm_3"}
};

class BipedSingleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    BipedSingleStepPeriodicityConstraints() = default;

    // Constructor
    BipedSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                          std::shared_ptr<BipedConstrainedInverseDynamics> dcidPtr_input,
                                          const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~BipedSingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<BipedConstrainedInverseDynamics> dcidPtr_;

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

}; // namespace Biped
}; // namespace RAPTOR

#endif // Biped_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
