#include "QRDecompositionSolver.h"

namespace RAPTOR {

QRDecompositionSolver::QRDecompositionSolver(const Model& model_input) {
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    phi = Eigen::VectorXd::Zero(10 * modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1;
        phi.segment<10>(10 * i) = 
            modelPtr_->inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }

    ObservationMatrix = MatX::Zero(0, 0);
}

void QRDecompositionSolver::generateRandomObservation(const int numInstances) {
    if (modelPtr_ == nullptr) {
        throw std::runtime_error("Model pointer is not initialized yet!");
    }
    
    ObservationMatrix = MatX::Zero(numInstances * modelPtr_->nv, phi.size());

    // Compute regressor matrix for random joint configurations
    for (int i = 0; i < numInstances; i++) {
        VecX q = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;
        VecX v = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;
        VecX a = 2 * M_PI * VecX::Random(modelPtr_->nv).array() - M_PI;

        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, 
            q, v, a);
            
        ObservationMatrix.middleRows(i * modelPtr_->nv, modelPtr_->nv) = 
            dataPtr_->jointTorqueRegressor;
    }
}

void QRDecompositionSolver::computeRegroupMatrix(const double eps) {
    if (ObservationMatrix.rows() == 0 ||
        ObservationMatrix.cols() == 0) {
        throw std::runtime_error("Observation matrix is not initialized yet!");
    }

    // Perform QR decomposition with column pivoting for stability
    Eigen::ColPivHouseholderQR<MatX> qr(ObservationMatrix);

    // Get permutation matrix
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();

    // Redefine the permutation matrix meaning, same as MATLAB
    // perm.indices()[i] = j old matrix ith col move to new matrix jth col
    // M[i] = j old jth col is M ith col
    std::vector<int> M(perm.indices().size());
    for (int i = 0; i < perm.indices().size(); i++) {
        M[perm.indices()[i]] = i;
    }

    // Get the rank of the matrix
    const int rankW = qr.rank();

    // Extract and sort independent cols
    std::vector<int> independent_idx(M.begin(), M.begin() + rankW);
    std::sort(independent_idx.begin(), independent_idx.end());

    // Extract and sort dependent cols
    std::vector<int> dependent_idx(M.begin() + rankW, M.end());
    std::sort(dependent_idx.begin(), dependent_idx.end());

    // Initialize Aid matrix
    Aid = MatX::Zero(ObservationMatrix.cols(), rankW);
    for (int i = 0; i < rankW; i++) {
        Aid(independent_idx[i], i) = 1.0;
    }

    // Initialize Ad matrix
    Ad = MatX::Zero(ObservationMatrix.cols(), ObservationMatrix.cols() - rankW);
    for (int i = 0; i < static_cast<int>(dependent_idx.size()); i++) {
        Ad(dependent_idx[i], i) = 1.0;
    }

    // Combine Aid and Ad to form the full selection matrix A
    A.resize(ObservationMatrix.cols(), ObservationMatrix.cols());
    A << Aid, Ad;

    // Perform QR decomposition on ObservationMatrix * A
    Eigen::ColPivHouseholderQR<MatX> qr_A(ObservationMatrix * A);
    MatX R_full = qr_A.matrixQR().triangularView<Eigen::Upper>();

    // Extract indep R1 and dep R2 
    const MatX& R1 = R_full.topLeftCorner(rankW, rankW);
    const MatX& R2 = R_full.topRightCorner(rankW, ObservationMatrix.cols() - rankW);

    // Compute Kd = R1^{-1} * R2 using a stable solver
    Kd = R1.colPivHouseholderQr().solve(R2);

    // Zero out elements in Kd that are below the threshold for numerical stability
    for (int i = 0; i < Kd.rows(); i++) {
        for (int j = 0; j < Kd.cols(); j++) {
            if (std::abs(Kd(i, j)) < eps) {
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
    MatX KG_(ObservationMatrix.cols(), ObservationMatrix.cols());
    KG_.setZero();
    KG_.topLeftCorner(rankW, rankW) = MatX::Identity(rankW, rankW);
    KG_.topRightCorner(rankW, ObservationMatrix.cols() - rankW) = -Kd;
    KG_.bottomRightCorner(ObservationMatrix.cols() - rankW, ObservationMatrix.cols() - rankW) = 
        MatX::Identity(ObservationMatrix.cols() - rankW, ObservationMatrix.cols() - rankW);

    Ginv = A * KG_;

    // Set the dimensions of independent and dependent parameters
    dim_id = phi_id.size();
    dim_d = phi_d.size();

    // Compute RegroupMatrix
    RegroupMatrix = Aid + Ad * Kd.transpose();
}

}; // namespace RAPTOR