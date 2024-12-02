#include "EndEffectorRegressorConditionNumber.h"

namespace RAPTOR {

EndEffectorRegressorConditionNumber::EndEffectorRegressorConditionNumber(
    std::shared_ptr<Trajectories>& trajPtr_input, 
    std::shared_ptr<RegressorInverseDynamics>& ridPtr_input) :
    trajPtr_(trajPtr_input),
    ridPtr_(ridPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void EndEffectorRegressorConditionNumber::compute(const VecX& z, 
                                                  bool compute_derivatives,
                                                  bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("EndEffectorRegressorConditionNumber does not support hessian computation");
    }

    // this will be called in ridPtr_->compute
    // trajPtr_->compute(z, compute_derivatives, compute_hessian);

    ridPtr_->compute(z, compute_derivatives, compute_hessian);

    // get the regressor corresponding to the end effector
    const int startCol = (ridPtr_->NB - 1) * 10; 
    const MatX& EndeffectorY = ridPtr_->Y.middleCols(startCol, 10);
  
    Eigen::JacobiSVD<MatX> svd(EndeffectorY, 
                               Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VecX& singularValues = svd.singularValues();
    const MatX& U = svd.matrixU();
    const MatX& V = svd.matrixV();

    const size_t lastRow = singularValues.size() - 1;
    const double& sigmaMax = singularValues(0);
    const double& sigmaMin = singularValues(lastRow);

    // log of the condition number in 2-norm 
    // (ratio between the largest and the smallest singular values)
    f = std::log(sigmaMax) - std::log(sigmaMin);

    if (compute_derivatives) {
        // refer to (17) in https://j-towns.github.io/papers/svd-derivative.pdf
        // for analytical gradient of singular values
        for (int i = 0; i < trajPtr_->varLength; i++) {
            const int startCol = (ridPtr_->NB - 1) * 10; 
            const MatX& gradEndeffectorY = ridPtr_->pY_pz(i).middleCols(startCol, 10);
            
            const double gradSigmaMax = U.col(0).transpose()       * gradEndeffectorY * V.col(0);
            const double gradSigmaMin = U.col(lastRow).transpose() * gradEndeffectorY * V.col(lastRow);

            grad_f(i) = gradSigmaMax / sigmaMax - gradSigmaMin / sigmaMin;
        }
    }
}

}; // namespace RAPTOR