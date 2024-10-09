#include "RegroupedLMIConstraints.h"

namespace RAPTOR {

RegroupedLMIConstraints::RegroupedLMIConstraints(const std::shared_ptr<QRDecompositionSolver>& qrSolverPtr_input,
                                                 const int num_links_input,
                                                 const int varLength) :
    qrSolverPtr_(qrSolverPtr_input),
    num_links(num_links_input) {
    if (varLength < qrSolverPtr_->dim_id + qrSolverPtr_->dim_d) {
        std::cerr << varLength << std::endl;
        std::cerr << qrSolverPtr_->dim_id + qrSolverPtr_->dim_d << std::endl;
        throw std::invalid_argument("Error: variable length is too short");
    }

    if (num_links * 10 != qrSolverPtr_->dim_id + qrSolverPtr_->dim_d) {
        std::cerr << num_links << std::endl;
        std::cerr << qrSolverPtr_->dim_id + qrSolverPtr_->dim_d << std::endl;
        throw std::invalid_argument("Error: number of links does not match the regrouped inertial parameters");
    }

    lmiConstraintsPtr_ = std::make_shared<LMIConstraints>(num_links, varLength);

    m = lmiConstraintsPtr_->m;
    initialize_memory(m, varLength);
}

void RegroupedLMIConstraints::compute(const VecX& z, 
                                      bool compute_derivatives,
                                      bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("Error: Hessian computation is not supported for LMI constraints");
    }        

    const int dim = qrSolverPtr_->dim_id + qrSolverPtr_->dim_d;

    // assume the inertial parameters start from the beginning of z
    // recover the original inertial parameters
    const VecX& phi = z.head(dim);
    VecX phi_recovered = qrSolverPtr_->Ginv * phi;

    // compute the LMI constraints
    lmiConstraintsPtr_->compute(phi_recovered, 
                                compute_derivatives, 
                                compute_hessian);

    g = lmiConstraintsPtr_->g;

    if (compute_derivatives) {
        pg_pz.setZero();
        pg_pz.leftCols(dim) =
            lmiConstraintsPtr_->pg_pz.leftCols(dim) * qrSolverPtr_->Ginv;
    }
}

void RegroupedLMIConstraints::compute_bounds() {
    g_lb = VecX::Zero(m);
    g_ub = VecX::Constant(m, 1e19);
}

void RegroupedLMIConstraints::print_violation_info() {
    for (int i = 0; i < num_links; i++) {
        for (int j = 0; j < 4; j++) {
            if (g(i * 4 + j) < 0) {
                std::cout << "RegroupedLMIConstraints.cpp: determinant of submatrix " 
                          << j 
                          << " for link " 
                          << i 
                          << " is violated: " 
                          << g(i * 4 + j) 
                          << " < 0" 
                          << std::endl;
            }
        }
    }
}

}; // namespace RAPTOR