#ifndef JOINTLIMITS_H
#define JOINTLIMITS_H

#include <memory>

#include "Constraints.h"
#include "Trajectories.h"

namespace IDTO {

class JointLimits : public Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using SpaMat = Eigen::SparseMatrix<double>;

    // Constructor
    JointLimits() = default;

    // Constructor
    JointLimits(std::unique_ptr<Trajectories> trajPtr_input, 
                const VecX& lowerLimits_input, 
                const VecX& upperLimits_input);

    // Destructor
    ~JointLimits() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bound
    virtual void compute_lb() override;

        // compute constraints upper bound
    virtual void compute_ub() override;

    // class variables:
    std::shared_ptr<Trajectories> trajPtr_;

    VecX lowerLimits;
    VecX upperLimits;

        // compute results are stored here
    VecX g;
    SpaMat pg_pz;

    VecX g_lb;
    VecX g_ub;
};

}; // namespace IDTO

#endif // JOINTLIMITS_H
