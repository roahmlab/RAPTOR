
#ifndef TALOS_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define TALOS_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

#include "Constraints.h"
#include "TalosConstrainedInverseDynamics.h"
#include "RectangleSurfaceContactConstraints.h"

#include <memory>
#include <string>

namespace RAPTOR {
namespace Talos {    

// ... ==> ... 
const std::string JOINT_MAP[NUM_INDEPENDENT_JOINTS][2] = {
    {"leg_left_1_joint", "leg_right_1_joint"},
    {"leg_left_2_joint", "leg_right_2_joint"},
    {"leg_left_3_joint", "leg_right_3_joint"},
    {"leg_left_4_joint", "leg_right_4_joint"},
    {"leg_left_5_joint", "leg_right_5_joint"},
    {"leg_left_6_joint", "leg_right_6_joint"},
    {"leg_right_1_joint", "leg_left_1_joint"},
    {"leg_right_2_joint", "leg_left_2_joint"},
    {"leg_right_3_joint", "leg_left_3_joint"},
    {"leg_right_4_joint", "leg_left_4_joint"},
    {"leg_right_5_joint", "leg_left_5_joint"},
    {"leg_right_6_joint", "leg_left_6_joint"}
};

class TalosSingleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::ModelTpl<double>;
    using Data = pinocchio::DataTpl<double>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TalosSingleStepPeriodicityConstraints() = default;

    // Constructor
    TalosSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                          std::shared_ptr<TalosConstrainedInverseDynamics> dcidPtr_input,
                                          const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~TalosSingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<TalosConstrainedInverseDynamics> dcidPtr_;

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

}; // namespace Talos
}; // namespace RAPTOR

#endif // TALOS_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
