#ifndef JOINTLIMITS_H
#define JOINTLIMITS_H

#include <memory>

#include "Constraints.h"
#include "Trajectories.h"
#include "Utils.h"

namespace IDTO {

class JointLimits : public Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    JointLimits() = default;

    // Constructor
    JointLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                const VecX& lowerLimits_input, 
                const VecX& upperLimits_input,
                const bool wrapToPiOrNot_input = false);

    // Destructor
    ~JointLimits() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override;

        // print violation information
    virtual void print_violation_info() override;

    // class variables:
    std::shared_ptr<Trajectories> trajPtr_;

    VecX lowerLimits;
    VecX upperLimits;
    
    bool wrapToPiOrNot = false;
};

}; // namespace IDTO

#endif // JOINTLIMITS_H
