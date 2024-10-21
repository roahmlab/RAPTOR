#ifndef LMI_CONSTRAINTS_H
#define LMI_CONSTRAINTS_H

#include "Constraints.h"
#include "pinocchio/spatial/symmetric3.hpp"

namespace RAPTOR {

// This is the base (abstract) class for all constraints
class LMIConstraints: public Constraints {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Mat4 = Eigen::Matrix4d;
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Symmetric3 = pinocchio::Symmetric3Tpl<double>;

    // Constructor
    LMIConstraints() = default;

    LMIConstraints(const int num_links_input,
                   const int varLength);

    // Destructor
    ~LMIConstraints() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    virtual void compute_bounds() override;

    virtual void print_violation_info() override;

    // class members:
    int num_links = 0;

    Eigen::Array<MatX, 1, 10> dLMIdz;
};

}; // namespace RAPTOR

#endif // LMI_CONSTRAINTS_H
