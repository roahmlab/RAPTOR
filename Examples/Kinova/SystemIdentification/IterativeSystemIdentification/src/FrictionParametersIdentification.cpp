#include "FrictionParametersIdentification.h"

namespace RAPTOR {

// [TNLP_set_parameters]
bool FrictionParametersIdentification::set_parameters(
    const Model& model_input,
    const std::shared_ptr<MatX>& posPtr_input,
    const std::shared_ptr<MatX>& velPtr_input,
    const std::shared_ptr<MatX>& accPtr_input,
    const std::shared_ptr<MatX>& torquePtr_input,
    const bool include_offset_input
) 
{
    enable_hessian = true;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    Nact = modelPtr_->nv;

    posPtr_ = posPtr_input;
    velPtr_ = velPtr_input;
    accPtr_ = accPtr_input;
    torquePtr_ = torquePtr_input;

    if (posPtr_->rows() != velPtr_->rows() || 
        posPtr_->rows() != accPtr_->rows() || 
        posPtr_->rows() != torquePtr_->rows() ||
        posPtr_->rows() != modelPtr_->nv) {
        throw std::invalid_argument("Error: input data matrices have different number of rows");
    }

    if (posPtr_->cols() != velPtr_->cols() || 
        posPtr_->cols() != accPtr_->cols() || 
        posPtr_->cols() != torquePtr_->cols()) {
        throw std::invalid_argument("Error: input data matrices have different number of columns");
    }

    N = posPtr_->cols();

    include_offset = include_offset_input;

    // compute nominal torque from data
    tau_inertials = MatX::Zero(Nact, N);
    for (Index i = 0; i < N; i++) {
        const VecX& q = posPtr_->col(i);
        const VecX& v = velPtr_->col(i);
        const VecX& a = accPtr_->col(i);

        pinocchio::rnea(*modelPtr_, *dataPtr_, q, v, a);

        tau_inertials.col(i) = dataPtr_->tau;
    } 

    // Random intialized the start point
    x0 = 10 * Eigen::VectorXd::Random(Nact * (include_offset ? 4 : 3)).array() + 10;

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
// returns some info about the nlp
bool FrictionParametersIdentification::get_nlp_info(
    Index& n, 
    Index& m,
    Index& nnz_jac_g, 
    Index& nnz_h_lag,
    IndexStyleEnum& index_style
)
{
    // number of decision variables
    n = Nact * (include_offset ? 4 : 3);
    numVars= n;

    // number of constraints
    m = 0;
    numCons= m;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // Use C-style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool FrictionParametersIdentification::get_bounds_info(
    Index n, 
    Number* x_l, 
    Number* x_u,
    Index m, 
    Number* g_l, 
    Number* g_u
)
{
    // static friction >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[i] = 0.0;
        x_u[i] = 50.0;
    }

    // damping >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[Nact + i] = 0.0;
        x_u[Nact + i] = 50.0;
    }

    // armature >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[2 * Nact + i] = 0.0;
        x_u[2 * Nact + i] = 50.0;
    }

    // static friction offset
    if (include_offset) {
        for (Index i = 0; i < Nact; i++) {
            x_l[3 * Nact + i] = -50.0;
            x_u[3 * Nact + i] = 50.0;
        }
    }

    // no bounds on constraints since m = 0
    // for (Index i = 0; i < m; i++) {
    //     g_l[i] = -1e19;
    //     g_u[i] = 0;
    // }

    return true;
}
// [get_bounds_info]

// [TNLP_eval_f]
// returns the value of the objective function
bool FrictionParametersIdentification::eval_f(
    Index n, 
    const Number* x, 
    bool new_x, 
    Number& obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX friction = z.head(Nact);
    VecX damping = z.segment(Nact, Nact);
    VecX armature = z.segment(2 * Nact, Nact);
    VecX offset = VecX::Zero(Nact);
    if (include_offset) {
        offset = z.tail(Nact);
    }
    
    obj_value = 0;

    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);
        const VecX& tau = torquePtr_->col(i);
        const VecX& tau_inertial = tau_inertials.col(i);

        VecX total_friction_force = 
            friction.cwiseProduct(q_d.cwiseSign()) +
            damping.cwiseProduct(q_d) +
            armature.cwiseProduct(q_dd) +
            offset;

        VecX tau_estimated = tau_inertial + total_friction_force;

        obj_value += 0.5 * (tau_estimated - tau).squaredNorm();
    }

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool FrictionParametersIdentification::eval_grad_f(
    Index n, 
    const Number* x, 
    bool new_x, 
    Number* grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);
    
    VecX friction = z.head(Nact);
    VecX damping = z.segment(Nact, Nact);
    VecX armature = z.segment(2 * Nact, Nact);
    VecX offset = VecX::Zero(Nact);
    if (include_offset) {
        offset = z.tail(Nact);
    }
    
    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);
        const VecX& tau = torquePtr_->col(i);
        const VecX& tau_inertial = tau_inertials.col(i);

        VecX total_friction_force = 
            friction.cwiseProduct(q_d.cwiseSign()) +
            damping.cwiseProduct(q_d) +
            armature.cwiseProduct(q_dd) +
            offset;

        VecX tau_estimated = tau_inertial + total_friction_force;
        VecX tau_diff = tau_estimated - tau;

        grad_f_vec.head(Nact) += tau_diff.cwiseProduct(q_d.cwiseSign());
        grad_f_vec.segment(Nact, Nact) += tau_diff.cwiseProduct(q_d);
        grad_f_vec.segment(2 * Nact, Nact) += tau_diff.cwiseProduct(q_dd);
        if (include_offset) {
            grad_f_vec.tail(Nact) += tau_diff;
        }
    }           

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i);
    }

    return true;
}

// [TNLP_eval_hess_f]
// return the hessian of the objective function hess_{x} f(x) as a dense matrix
bool FrictionParametersIdentification::eval_hess_f(
   Index         n,
   const Number* x,
   bool          new_x,
   MatX&         hess_f
)
{
    if (n != numVars) {
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    hess_f = MatX::Zero(n, n);

    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);

        hess_f.diagonal().head(Nact) += q_d.cwiseSign().cwiseAbs2();
        hess_f.block(0, Nact, Nact, Nact).diagonal().array() += q_d.cwiseSign().array() * q_d.array();
        hess_f.block(0, 2 * Nact, Nact, Nact).diagonal().array() += q_d.cwiseSign().array() * q_dd.array();
        if (include_offset) {
            hess_f.block(0, 3 * Nact, Nact, Nact).diagonal() += q_d.cwiseSign();
        }

        hess_f.diagonal().segment(Nact, Nact) += q_d.cwiseAbs2();
        hess_f.block(Nact, 2 * Nact, Nact, Nact).diagonal().array() += q_d.array() * q_dd.array();
        if (include_offset) {
            hess_f.block(Nact, 3 * Nact, Nact, Nact).diagonal() += q_d;
        }

        hess_f.diagonal().segment(2 * Nact, Nact) += q_dd.cwiseAbs2();
        if (include_offset) {
            hess_f.block(2 * Nact, 3 * Nact, Nact, Nact).diagonal() += q_dd;
        }

        if (include_offset) {
            hess_f.diagonal().tail(Nact) += VecX::Ones(Nact);
        }
    }

    return true;
}
// [TNLP_eval_hess_f]

}; // namespace RAPTOR