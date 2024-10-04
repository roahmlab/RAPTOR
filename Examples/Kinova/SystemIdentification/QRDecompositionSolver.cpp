#include "QRDecompositionSolver.h"

namespace RAPTOR {

// Default constructor
QRDecompositionSolver::QRDecompositionSolver() : eps(1e-8), dim_id(0), dim_d(0) {}

// Parameterized constructor
QRDecompositionSolver::QRDecompositionSolver(const MatX& ObservationMatrix, const VecX& InputParams)
    : eps(1e-8), dim_id(0), dim_d(0) {
    compute(ObservationMatrix, InputParams);
}

// Destructor
QRDecompositionSolver::~QRDecompositionSolver() {}

// load data (robot model and generate ObservationMatrix and InputParams)
void QRDecompositionSolver::getData() {
    // Load robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);


    Eigen::VectorXd X0_1_ = Eigen::VectorXd::Zero(10*model.nv);
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1;
        // intertia 
         X0_1_.segment<10>(10 * i) = model.inertias[pinocchio_joint_id].toDynamicParameters();
    }

    InputParams= X0_1_;

    // Create a trajectory
    int N = 1000;  // Number of time steps
    double T = 1000.0;  // Total time
    int degree = 5;  // Degree of the polynomial

    // Assuming Polynomials constructor 
    trajPtr_ = std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Assuming RegressorInverseDynamics constructor
    ridPtr_ = std::make_shared<RegressorInverseDynamics>(model, trajPtr_);

    // Generate random joint positions, velocities, and accelerations (not accurate)
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    VecX z = M_2_PI * VecX::Random(trajPtr_->varLength).array() - M_PI;

    // Assuming compute is a member function of RegressorInverseDynamics
    ridPtr_->compute(z, false);

    // Retrieve ObservationMatrix and InputParams
    MatX ObservationMatrix = ridPtr_->Y.leftCols(model.nv *10);

    // InputParams = VecX::Zero(91); 
    // InputParams = VecX::Ones(70); 

    // Perform QR Decomposition
    compute(ObservationMatrix, InputParams);
}

// Compute QR Decomposition-based regrouping
void QRDecompositionSolver::compute(const MatX& W, const VecX& InputParams) {
    // Perform QR decomposition with column pivoting for stability
    Eigen::ColPivHouseholderQR<MatX> qr(W);

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

    int p = W.cols();

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

    // Perform QR decomposition on W * A
    Eigen::HouseholderQR<MatX> qr_A(W * A);
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
    pi_id = Aid.transpose() * InputParams;
    pi_d = Ad.transpose() * InputParams;

    // Compute base inertial parameters
    Beta = pi_id + Kd * pi_d;

    // Compute Ginv
    MatX KG_(p, p);
    KG_.setZero();
    KG_.topLeftCorner(rankW, rankW) = MatX::Identity(rankW, rankW);
    KG_.topRightCorner(rankW, p - rankW) = -Kd;
    KG_.bottomRightCorner(p - rankW, p - rankW) = MatX::Identity(p - rankW, p - rankW);
    Ginv = A * KG_;

    // Set the dimensions of independent and dependent parameters
    dim_id = pi_id.size();
    dim_d = pi_d.size();

    // Compute RegroupMatrix
    RegroupMatrix = Aid + Ad * Kd.transpose();
}

}; // namespace RAPTOR
