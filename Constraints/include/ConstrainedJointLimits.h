#ifndef CONSTRAINEDJOINTLIMITS_H
#define CONSTRAINEDJOINTLIMITS_H

#include "JointLimits.h"
#include "DynamicsConstraints.h"

namespace IDTO {

class ConstrainedJointLimits : public JointLimits {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using SpaMatX = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    // Constructor
    ConstrainedJointLimits() = default;

    // Constructor
    ConstrainedJointLimits(std::unique_ptr<Trajectories> trajPtr_input, 
                           std::unique_ptr<DynamicsConstraints> dcPtr_input, 
                           const VecX& lowerLimits_input, 
                           const VecX& upperLimits_input);

    // Destructor
    ~ConstrainedJointLimits() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bound
    void compute_lb() override;

        // compute constraints upper bound
    void compute_ub() override;

    // class variables:
    int NB = 0;
    std::shared_ptr<DynamicsConstraints> dcPtr_;
};

}; // namespace IDTO

#endif // CONSTRAINEDJOINTLIMITS_H
