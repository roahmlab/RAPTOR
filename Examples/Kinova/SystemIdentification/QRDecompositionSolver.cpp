#include "QRDecompositionSolver.h"
#include <algorithm>
#include <limits>
#include <ctime>
#include <cstdlib>
#include <cmath>

namespace RAPTOR {
namespace Kinova {

QRDecompositionSolver::QRDecompositionSolver(){
    getdata();
}

QRDecompositionSolver::QRDecompositionSolver(const MatX& ObservationMatrix, const VecX& InParam) {
    compute(ObservationMatrix, InParam);
}

QRDecompositionSolver::~QRDecompositionSolver() {}


void QRDecompositionSolver::getdata(){
    // Load robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // Create a trajectory
    int N = 1000;  // number of time steps
    double T = 1000.0;  // total time
    int degree = 5;  // degree of the polynomial
    trajPtr_ = std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);
    ridPtr_ = std::shared_ptr<RegressorInverseDynamics>(model, trajPtr_);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    VecX z = 2 * M_PI * VecX::Random(trajPtr->varLength).array() - M_PI;

    ridPtr_.compute(z, false);
    MatX ObservationMatrix = ridPtr_.Y;
    VecX InParam = VecX::Zero(91);

    // compute
    compute(ObservationMatrix, InParam);
}

void RDecompositionSolver::compute(const MatX& W, const VecX& InParam){
    

    // QR decomposition with column pivoting, improve stability
    Eigen::ColPivHouseholderQR<MatX> qr(W);

    //  Permutation matrix
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();

    // Redefine the permutation matrix meaning, same as matlab
    // perm.indices()[i] = j   old matrix ith col move to new matrix jth col
    // M[i] = j  old  j jth col is M ith col
    std::vector<int> M(perm.indices().size());
    for (int i =0; i< perm.indices().size(); ++i){
        M[perm.indices()[i]]= i;
    }


    int rankW = qr.rank();
    std::vector<int> idx(M.begin(), M.begin() + rankW);
    std::sort(idx.begin(), idx.end());
    std::vector<int> idx_(M.begin() + rankW, M.end());
    std::sort(idx_.begin(), idx_.end());

    int p = W.cols();
    Aid = MatX::Zero(p, rankW);
    for(int i = 0; i < rankW; ++i) {
        Aid(idx[i], i) = 1.0;
    }

    Ad = MatX::Zero(p, p - rankW);
    for(int i = 0; i < static_cast<int>(idx_.size()); ++i) {
        Ad(idx_[i], i) = 1.0;
    }

    // Combine Aid and Ad 
    MatX A(p,p);
    A << Aid, Ad;

    // QR decomposition on W * A
    Eigen::HouseholderQR<MatX> qr_A(W * A);
    MatX R_full = qr_A.matrixQR().triangularView<Eigen::Upper>();

    // Extract indep R1 (upper-left block) and dep R2 (upper-right block)
    MatX R1 = R_full.topLeftCorner(rankW, rankW);
    MatX R2 = R_full.topRightCorner(rankW, p - rankW);

    // kd = (R1^-1)*R2  b*(p-b)
    Kd = R1.solve(R2) 

    // get rid of extremely small numbers for more numerical stability
    double threshold = std::sqrt(eps);
    for(int i = 0; i < Kd.rows(); ++i) {
        for(int j = 0; j < Kd.cols(); ++j) {
            if(std::abs(Kd(i, j)) < threshold) {
                Kd(i, j) = 0.0;
            }
        }
    }

    // Get independant and dependant parameters
    VecX pi_id = Aid.transpose() * InParam;
    VecX pi_d = Ad.transpose() * InParam;

    // Compute base inertial parameters
    beta = pi_id + Kd * pi_d;

    // compute Ginv
    MatX KG_(p,p);
    KG_.setZero();
    KG_.topLeftCorner(rankW, rankW) = MatX::Identity(rankW, rankW);
    KG_.topRightCorner(rankW, p - rankW) = -Kd;
    KG_.bottomRightCorner(p - rankW, p - rankW) = MatX::Identity(p - rankW, p - rankW);
    Ginv = A * KG_;

    // Set the dimensions of independent and dependent parameters
    dim_id = pi_id.size();
    dim_d = pi_d.size();

    // Compute RegroupMatrix
    RegroupMatrix = Aid + Ad * Kd.transpose()
}

}; // namespace Kinova
}; // namespace RAPTOR
