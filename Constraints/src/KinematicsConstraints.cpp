#include "KinematicsConstraints.h"

namespace IDTO {

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                       const Eigen::VectorXi& jtype_input) : 
    trajPtr_(trajPtr_input),
    jtype(jtype_input) {
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    
}

void KinematicsConstraints::compute(const VecX& z, bool compute_derivatives) {

}

}; // namespace IDTO