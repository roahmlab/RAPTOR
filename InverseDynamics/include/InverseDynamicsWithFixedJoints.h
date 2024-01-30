#ifndef INVERSEDYNAMICS_WITHFIXEDJOINTS_H
#define INVERSEDYNAMICS_WITHFIXEDJOINTS_H

#include "InverseDynamics.h"

namespace IDTO {

#define FINITE_DIFFERENCE_DELTA 1e-8

class InverseDynamicsWithFixedJoints : public InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Force = Data::Force;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    InverseDynamicsWithFixedJoints() = default;

    // Constructor
    InverseDynamicsWithFixedJoints(const Model& model_input, 
                                   std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~InverseDynamicsWithFixedJoints() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true) override;   

    // class members:
        // the contact wrench at the last joint, which is assumed to be fixed
    Eigen::Array<Force, 1, Eigen::Dynamic> lambda;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pz;
    
};

}; // namespace IDTO

#endif // INVERSEDYNAMICS_WITHFIXEDJOINTS_H
