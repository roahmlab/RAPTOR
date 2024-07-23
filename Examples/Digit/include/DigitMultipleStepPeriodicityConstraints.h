
#ifndef DIGIT_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H
#define DIGIT_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "DigitSingleStepOptimizer.h"

namespace RAPTOR {
namespace Digit {

class DigitMultipleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitMultipleStepPeriodicityConstraints() = default;

    // Constructor
    DigitMultipleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& currTrajPtr_input,
                                            std::shared_ptr<Trajectories>& nextTrajPtr_input,
                                            std::shared_ptr<DigitConstrainedInverseDynamics> currDcidPtr_input,
                                            std::shared_ptr<DigitConstrainedInverseDynamics> nextDcidPtr_input,
                                            const frictionParams& fp_input);

    // Destructor
    ~DigitMultipleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& currTrajPtr_;
    std::shared_ptr<Trajectories>& nextTrajPtr_;
    std::shared_ptr<DigitConstrainedInverseDynamics> currDcidPtr_;
    std::shared_ptr<DigitConstrainedInverseDynamics> nextDcidPtr_;

    frictionParams fp;

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

}; // namespace Digit
}; // namespace RAPTOR

#endif // DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
