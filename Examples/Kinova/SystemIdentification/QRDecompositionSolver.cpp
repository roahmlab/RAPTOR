#include "QRDecompositionSolver.h"

namespace RAPTOR {

QRDecompositionSolver::QRDecompositionSolver(const Model& model) {
    modelPtr_ = std::make_shared<Model>(model);
    dataPtr_ = std::make_shared<Data>(model);

    phi = Eigen::VectorXd::Zero(10 * modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1;
        phi.segment<10>(10 * i) = modelPtr_->inertias[pinocchio_joint_id].toDynamicParameters();
    }

    ObservationMatrix = MatX::Zero(0, 0);
}

void QRDecompositionSolver::generateRandomObservation(const int numInstances) {
    if (modelPtr_ == nullptr) {
        throw std::runtime_error("Model pointer is not initialized yet!");
    }
    
    ObservationMatrix = MatX::Zero(numInstances * modelPtr_->nv, phi.size());

    for (int i = 0; i < numInstances; i++) {
        // Generate random joint configuration
        VecX q = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;
        VecX v = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;
        VecX a = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;

        // Compute regressor
        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, 
            q, v, a);

        // Fill in the observation matrix
        ObservationMatrix.middleRows(i * modelPtr_->nv, modelPtr_->nv) = 
            dataPtr_->jointTorqueRegressor;
    }
}

void QRDecompositionSolver::computeRegroupMatrix() {
    if (ObservationMatrix.rows() == 0 ||
        ObservationMatrix.cols() == 0) {
        throw std::runtime_error("Observation matrix is not initialized yet!");
    }

    // Perform QR decomposition with column pivoting for stability
    Eigen::ColPivHouseholderQR<MatX> qr(ObservationMatrix);

    // Get permutation matrix
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();

    // Redefine the permutation matrix meaning, same as MATLAB
    // perm.indices()[i] = j   old matrix ith col move to new matrix jth col
    // M[i] = j  old jth col is M ith col
    std::vector<int> M(perm.indices().size());
    for (int i = 0; i < perm.indices().size(); ++i) {
        M[perm.indices()[i]] =i;
    }

    // Get the rank of the matrix
    int rankW = qr.rank();

    // Extract and sort independent cols
    std::vector<int> idx(M.begin(), M.begin() + rankW);
    std::sort(idx.begin(), idx.end());

    // Extract and sort dependent cols
    std::vector<int> idx_(M.begin() + rankW, M.end());
    std::sort(idx_.begin(), idx_.end());

    int p = ObservationMatrix.cols();

    // Initialize Aid matrix
    Aid = MatX::Zero(p, rankW);
    for (int i = 0; i < rankW; ++i) {
        Aid(idx[i], i) = 1.0;
    }

    // Initialize Ad matrix
    Ad = MatX::Zero(p, p - rankW);
    for (int i = 0; i < static_cast<int>(idx_.size()); ++i) {
        Ad(idx_[i], i) = 1.0;
    }

    // Combine Aid and Ad to form the full selection matrix A
    A.resize(p, p);
    A << Aid, Ad;

    // Perform QR decomposition on ObservationMatrix * A
    Eigen::HouseholderQR<MatX> qr_A(ObservationMatrix * A);
    MatX R_full = qr_A.matrixQR().triangularView<Eigen::Upper>();

    // Extract indep R1 and dep R2 
    MatX R1 = R_full.topLeftCorner(rankW, rankW);
    MatX R2 = R_full.topRightCorner(rankW, p - rankW);

    // Compute Kd = R1^{-1} * R2 using a stable solver
    Kd = R1.colPivHouseholderQr().solve(R2);

    // Zero out elements in Kd that are below the threshold for numerical stability
    double threshold = std::sqrt(eps);
    for (int i = 0; i < Kd.rows(); ++i) {
        for (int j = 0; j < Kd.cols(); ++j) {
            if (std::abs(Kd(i, j)) < threshold) {
                Kd(i, j) = 0.0;
            }
        }
    }

    // Get independent and dependent parameters
    phi_id = Aid.transpose() * phi;
    phi_d = Ad.transpose() * phi;

    // Compute base inertial parameters
    beta = phi_id + Kd * phi_d;

    // Compute Ginv
    MatX KG_(p, p);
    KG_.setZero();
    KG_.topLeftCorner(rankW, rankW) = MatX::Identity(rankW, rankW);
    KG_.topRightCorner(rankW, p - rankW) = -Kd;
    KG_.bottomRightCorner(p - rankW, p - rankW) = MatX::Identity(p - rankW, p - rankW);
    Ginv = A * KG_;

    // Set the dimensions of independent and dependent parameters
    dim_id = phi_id.size();
    dim_d = phi_d.size();

    // Compute RegroupMatrix
    RegroupMatrix = Aid + Ad * Kd.transpose();
}

}; // namespace RAPTOR
