#ifndef CONSTRAINEDJOINTLIMITS_H
#define CONSTRAINEDJOINTLIMITS_H

#include "JointLimits.h"
#include "DynamicsConstraints.h"

namespace RAPTOR {

class ConstrainedJointLimits : public JointLimits {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    ConstrainedJointLimits() = default;

    // Constructor
    ConstrainedJointLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                           std::shared_ptr<DynamicsConstraints>& dcPtr_input, 
                           const VecX& lowerLimits_input, 
                           const VecX& upperLimits_input);

    // Destructor
    ~ConstrainedJointLimits() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override;

    // class variables:
    int NB = 0;
    std::shared_ptr<DynamicsConstraints> dcPtr_;
};

}; // namespace RAPTOR

#endif // CONSTRAINEDJOINTLIMITS_H
