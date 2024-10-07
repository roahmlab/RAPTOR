#ifndef REGROUPED_LMI_CONSTRAINTS_H
#define REGROUPED_LMI_CONSTRAINTS_H

#include "LMIConstraints.h"
#include "QRDecompositionSolver.h"

namespace RAPTOR {

// This is the base (abstract) class for all constraints
class RegroupedLMIConstraints: public Constraints {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Mat4 = Eigen::Matrix4d;
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Symmetric3 = pinocchio::Symmetric3Tpl<double>;

    // Constructor
    RegroupedLMIConstraints() = default;

    RegroupedLMIConstraints(const std::shared_ptr<QRDecompositionSolver>& qrSolverPtr_input,
                            const int num_links_input,
                            const int varLength);

    // Destructor
    ~RegroupedLMIConstraints() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    virtual void compute_bounds() override;

    virtual void print_violation_info() override;

    // class members:
    int num_links = 0;

    std::shared_ptr<QRDecompositionSolver> qrSolverPtr_ = nullptr;
    std::shared_ptr<LMIConstraints> lmiConstraintsPtr_ = nullptr;
};

}; // namespace RAPTOR

#endif // REGROUPED_LMI_CONSTRAINTS_H
