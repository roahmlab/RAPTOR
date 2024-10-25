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
    // step 1: Perform QR decomposition with column pivoting for stability
    Eigen::ColPivHouseholderQR<MatX> qr(ObservationMatrix);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();
    const int rankW = qr.rank();
    A.resize(ObservationMatrix.cols(), ObservationMatrix.cols());
    A = perm;
    Aid = A.leftCols(rankW);
    Ad = A.rightCols(ObservationMatrix.cols()-rankW);

    // step 2: Extract indep R1 and dep R2
    MatX R_full = qr.matrixQR().triangularView<Eigen::Upper>();
    const MatX& R1 = R_full.topLeftCorner(rankW, rankW);
    const MatX& R2 = R_full.topRightCorner(rankW, ObservationMatrix.cols() - rankW);
    // sanity check on rank of R1
    if (R1.colPivHouseholderQr().rank() != rankW) {
        std::cerr << R1.colPivHouseholderQr().rank() << std::endl;
        std::cerr << rankW << std::endl;
        throw std::runtime_error("Rank of R1 is not full rank!");
    }

    // step 3: calculate the necessary varialbes 
    // Compute Kd = R1^{-1} * R2 using a stable solver
    Kd = R1.completeOrthogonalDecomposition().pseudoInverse() * R2;
    // Zero out elements in Kd that are below the threshold for numerical stability
    for (int i = 0; i < Kd.rows(); i++) {
        for (int j = 0; j < Kd.cols(); j++) {
            if (std::abs(Kd(i, j)) < eps) {
                Kd(i, j) = 0.0;
            }
        }
    }

    // Get independent and dependent and base inertial parameters
    phi_id = Aid.transpose() * phi;
    phi_d = Ad.transpose() * phi;
    dim_id = phi_id.size();
    dim_d = phi_d.size();
    beta = phi_id + Kd * phi_d;

    // Compute Ginv
    MatX KG_(ObservationMatrix.cols(), ObservationMatrix.cols());
    KG_.setZero();
    KG_.topLeftCorner(rankW, rankW) = MatX::Identity(rankW, rankW);
    KG_.topRightCorner(rankW, ObservationMatrix.cols() - rankW) = -Kd;
    KG_.bottomRightCorner(ObservationMatrix.cols() - rankW, ObservationMatrix.cols() - rankW) =
        MatX::Identity(ObservationMatrix.cols() - rankW, ObservationMatrix.cols() - rankW);
    Ginv = A * KG_;

    // Compute RegroupMatrix
    RegroupMatrix = Aid + Ad * Kd.transpose();
}
}; // namespace RAPTOR