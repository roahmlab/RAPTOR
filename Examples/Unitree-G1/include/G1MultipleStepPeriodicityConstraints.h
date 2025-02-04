
#ifndef G1_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H
#define G1_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "G1SingleStepOptimizer.h"

namespace RAPTOR {
namespace G1 {

class G1MultipleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    G1MultipleStepPeriodicityConstraints() = default;

    // Constructor
    G1MultipleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& currTrajPtr_input,
                                         std::shared_ptr<Trajectories>& nextTrajPtr_input,
                                         std::shared_ptr<G1ConstrainedInverseDynamics> currDcidPtr_input,
                                         std::shared_ptr<G1ConstrainedInverseDynamics> nextDcidPtr_input,
                                         const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~G1MultipleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& currTrajPtr_;
    std::shared_ptr<Trajectories>& nextTrajPtr_;
    std::shared_ptr<G1ConstrainedInverseDynamics> currDcidPtr_;
    std::shared_ptr<G1ConstrainedInverseDynamics> nextDcidPtr_;

    rectangleContactSurfaceParams fp;

        // intermediate variables updated in compute()
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
    MatX pg3_pz2;
    MatX pg4_pz;
    MatX pg4_pz2;
    MatX pg5_pz;
};

}; // namespace G1
}; // namespace RAPTOR

#endif // G1_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
