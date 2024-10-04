#include "FrictionParametersIdentification.h"
#include <iostream>
#include <cmath>

namespace RAPTOR {
namespace Kinova {


bool FrictionParametersIdentification::set_parameters(
        VecX Xf,
        int nLinks,
        VecX& Fest,
        bool include_friction_offset,
        double N,
        std::shared_ptr<RegressorInverseDynamics>& RegressorID
) 

{
    Xf_=Xf;
    nLinks_= nLinks;
    Fest_=Fest;
    include_friction_offset_=include_friction_offset; 
    N_=N;
    ridPtr_= RegressorID;

    enable_hessian = false;
    return true;
}

bool FrictionParametersIdentification::get_nlp_info(Index& n, Index& m,
                                        Index& nnz_jac_g, Index& nnz_h_lag,
                                        IndexStyleEnum& index_style)
{
    // αj and phiF parameters
    int n_alpha = nLinks_;
    // Fvj, Fcj, B 
    int n_phiF = nLinks_ * (include_friction_offset_ ? 3 : 2);
    n = n_alpha + n_phiF;
    numVars= n;

    // Number of constraints (none in this problem)
    // m = nLinks_;
    m=0;
    numCons= m;

    nnz_jac_g = 0;
    nnz_h_lag = n * (n + 1) / 2;

    // Use C-style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool FrictionParametersIdentification::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                           Index m, Number* g_l, Number* g_u)
{
    // αj >= 0
    for (Index i = 0; i < nLinks_; ++i) {
        x_l[i] = 0.0;
        x_u[i] = 1e19;
    }

    // phiF parameters bounds
    Index idx = nLinks_;
    for (Index j = 0; j < nLinks_; ++j) {
        // Fcj >= 0
        x_l[idx] = 0.0;
        x_u[idx] = 1e19;
        idx++;

        // Fvj >= 0
        x_l[idx] = 0.0;
        x_u[idx] = 1e19;
        idx++;

        // No bounds on Bj
        if (include_friction_offset_) {
            // No bounds on Bj
            x_l[idx] = -1e19;
            x_u[idx] = 1e19;
            idx++;
        }

        // for (Index i = 0; i < m; ++i)
        // {
        //     g_l[i] = -1e19;
        //     g_u[i] = 0;
        // }

    }

    return true;
}

bool FrictionParametersIdentification::get_starting_point(Index n, bool init_x, Number* x,
                                              bool init_z, Number* z_L, Number* z_U,
                                              Index m, bool init_lambda,
                                              Number* lambda)
{
    // Initialize x
    init_x = true;
    init_z = false;
    init_lambda = false;

    // Initialize αj
    for (Index i = 0; i < nLinks_; ++i) {
        x[i] = Xf_[i]; // Initial guess for αj
    }

    // Initialize phiF parameters
    Index idx = nLinks_;
    for (Index j = 0; j < nLinks_; ++j) {
        x[idx] = Xf_[idx]; 
        idx++;
        x[idx] = Xf_[idx]; 
        idx++;
        if (include_friction_offset_) {
            x[idx] = Xf_[idx]; 
            idx++;
        }
    }
    return true;
}

bool FrictionParametersIdentification::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    // Extract αj phiF
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX alpha = z.head(nLinks_); 
    VecX phiF = z.tail(n - nLinks_);

    // Compute the objective function value
    obj_value = 0.0;

    for (Index i = 0; i < N_; ++i) {
        const VecX& q_d = ridPtr_->trajPtr_q_d(i);
        for (Index j = 0; j < nLinks_; ++j) {
            double sign_qd = (q_d(j) >= 0) ? 1.0 : -1.0;
            double abs_qd = std::abs(q_d(j));
            double abs_qd_alpha = std::pow(abs_qd, alpha(j));

            // Get phiF parameters for joint j
            Index phiF_offset = j * (include_friction_offset_ ? 3 : 2);
            double Fcj = phiF(phiF_offset);
            double Fvj = phiF(phiF_offset + 1);
            double Bj =0.0;
            if (include_friction_offset_) {
                Bj = phiF(phiF_offset + 2);
            }

            // Compute friction model fj
            double fj = (Fcj + Fvj * abs_qd_alpha) * sign_qd + Bj;

            // Compute obj_value
            double diff = fj - Fest_(i*nLinks_ + j);
            obj_value += diff * diff;
        }
    }

    return true;
}

// Method to compute gradient of the objective function
bool FrictionParametersIdentification::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    // Extract αj phiF
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX alpha = z.head(nLinks_); 
    VecX phiF = z.tail(n - nLinks_);

    // Initialize gradient to zero
    VecX grad = VecX::Zero(n);

    // Compute gradient
    for (Index i = 0; i < N_; ++i) {
        const VecX& q_d = ridPtr_->trajPtr_q_d(i);
        for (Index j = 0; j < nLinks_; ++j) {
            double sign_qd = (q_d(j) >= 0) ? 1.0 : -1.0;
            double abs_qd = std::abs(q_d(j));
            double abs_qd_alpha = std::pow(abs_qd, alpha(j));
            //avoid log(0)
            double log_abs_qd = (abs_qd > 1e-6) ? std::log(abs_qd) : 0.0;

            // Get phiF parameters for joint j
            Index phiF_offset = j * (include_friction_offset_ ? 3 : 2);
            double Fcj = phiF(phiF_offset);
            double Fvj = phiF(phiF_offset + 1);
            double Bj = 0.0;
            if (include_friction_offset_) {
                Bj = phiF(phiF_offset + 2);
            }

            // Compute friction model fj
            double fj = (Fcj + Fvj * abs_qd_alpha) * sign_qd + Bj;
            double diff = fj - Fest_(i*nLinks_ + j);

            // dalpha 
            double dfj_dalpha = Fvj * abs_qd_alpha * log_abs_qd * sign_qd; 
            grad(j) += 2.0 * diff * dfj_dalpha;
            // dFvj
            grad(nLinks_ + phiF_offset + 1) += 2.0 * diff * abs_qd_alpha * sign_qd;
            // dFcj
            grad(nLinks_ + phiF_offset) += 2.0 * diff * sign_qd;
            // Bj
            if (include_friction_offset_) {
                grad(nLinks_ + phiF_offset + 2) += 2.0 * diff;
            }
        }
    }
    // Copy gradient to output
    Eigen::Map<VecX>(grad_f, n) = grad;


    return true;
}


// bool FrictionParametersIdentification::eval_g(Index n, const Number* x, bool new_x,
//                                   Index m, Number* g)
// {


//     return true;
// }


}; // namespace Kinova
}; // namespace RAPTOR