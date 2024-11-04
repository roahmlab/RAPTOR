
#ifndef TALOS_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H
#define TALOS_MULTIPLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "TalosSingleStepOptimizer.h"

namespace RAPTOR {
namespace Talos {

class TalosMultipleStepPeriodicityConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::DataTpl<double>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TalosMultipleStepPeriodicityConstraints() = default;

    // Constructor
    TalosMultipleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& currTrajPtr_input,
                                            std::shared_ptr<Trajectories>& nextTrajPtr_input,
                                            std::shared_ptr<TalosConstrainedInverseDynamics> currDcidPtr_input,
                                            std::shared_ptr<TalosConstrainedInverseDynamics> nextDcidPtr_input,
                                            const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~TalosMultipleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    void compute_bounds() final override;

    void print_violation_info() final override;

    // class members
    std::shared_ptr<Trajectories>& currTrajPtr_;
    std::shared_ptr<Trajectories>& nextTrajPtr_;
    std::shared_ptr<TalosConstrainedInverseDynamics> currDcidPtr_;
    std::shared_ptr<TalosConstrainedInverseDynamics> nextDcidPtr_;

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

}; // namespace Talos
}; // namespace RAPTOR

#endif // TALOS_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
