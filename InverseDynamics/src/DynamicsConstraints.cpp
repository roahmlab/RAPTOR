#include "DynamicsConstraints.h"

namespace IDTO {

DynamicsConstraints::DynamicsConstraints(const Model& model_input, int numDependentJoints) {
    modelPtr_ = std::make_unique<Model>(model_input);

    c = VecX::Zero(modelPtr_->nv);
    J = MatX::Zero(numDependentJoints, modelPtr_->nv);
    Jx_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTx_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    Jxy_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
}

}; // namespace IDTO