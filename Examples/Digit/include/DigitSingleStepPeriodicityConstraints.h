
#ifndef DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
#define DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H

#include "pinocchio/algorithm/crba.hpp"

#include "Constraints.h"
#include "DigitConstrainedInverseDynamics.h"

#include <memory>

namespace IDTO {
namespace Digit {    

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
                                          std::shared_ptr<DigitConstrainedInverseDynamics> dcidPtr_input);

    // Destructor
    ~DigitSingleStepPeriodicityConstraints() = default;

    // class methods
    void compute(const VecX& z, bool compute_derivatives = true) override;

    void compute_bounds() override;

    // class members
    std::shared_ptr<Trajectories>& trajPtr_;
    std::shared_ptr<DigitConstrainedInverseDynamics> dcidPtr_;

};

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_SINGLE_STEP_PERIODICITY_CONSTRAINTS_H
