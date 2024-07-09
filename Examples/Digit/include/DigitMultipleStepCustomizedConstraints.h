#ifndef DIGIT_MULTIPLE_STEP_CUSTOMIZED_CONSTRAINTS_H
#define DIGIT_MULTIPLE_STEP_CUSTOMIZED_CONSTRAINTS_H

#include "DigitCustomizedConstraints.h"

namespace IDTO {
namespace Digit {

class DigitMultipleStepCustomizedConstraints : public DigitCustomizedConstraints {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitMultipleStepCustomizedConstraints() = default;

    // Constructor
    DigitMultipleStepCustomizedConstraints(const Model& model_input,
                                           const Eigen::VectorXi& jtype_input,
                                           std::shared_ptr<Trajectories>& trajPtr_input,
                                           std::shared_ptr<DigitDynamicsConstraints>& dcPtr_input,
                                           const GaitParameters& gp_input,
                                           const int N_input,
                                           const int NSteps_input);

    // Destructor
    ~DigitMultipleStepCustomizedConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() final override;

    // class variables:
    int N = 0;
    int NSteps = 0;
};

} // namespace Digit
} // namespace IDTO

#endif // DIGIT_MULTIPLE_STEP_CUSTOMIZED_CONSTRAINTS_H
