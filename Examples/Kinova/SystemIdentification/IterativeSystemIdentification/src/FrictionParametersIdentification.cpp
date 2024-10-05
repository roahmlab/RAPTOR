#include "FrictionParametersIdentification.h"

namespace RAPTOR {
namespace Kinova {

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
    
    return true;
}

bool FrictionParametersIdentification::get_nlp_info(
    Index& n, 
    Index& m,
    Index& nnz_jac_g, 
    Index& nnz_h_lag,
    IndexStyleEnum& index_style
)
{
    n = Nact * (include_offset ? 4 : 3);
    numVars= n;

    m = 0;
    numCons= m;

    nnz_jac_g = 0;
    nnz_h_lag = 0;

    // Use C-style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

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

bool FrictionParametersIdentification::get_starting_point(
    Index n, 
    bool init_x, 
    Number* x,
    bool init_z, 
    Number* z_L, 
    Number* z_U,
    Index m, 
    bool init_lambda,
    Number* lambda
)
{
    // simply use zero as starting point
    for (Index i = 0; i < n; i++) {
        x[i] = 0;
    }

    return true;
}

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

        obj_value += (tau_estimated - tau).squaredNorm();
    }

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

// Method to compute gradient of the objective function
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

        // TODO: compute gradient of the objective function
    }

    return true;
}

bool FrictionParametersIdentification::eval_hess_f(
   Index         n,
   const Number* x,
   bool          new_x,
   MatX&         hess_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    // TODO: according to eval_grad_f, verify this
    // theoretically, this is a quadratic problem so the hessian should be idenity matrix
    hess_f = MatX::Identity(n, n);

    return true;
}

}; // namespace Kinova
}; // namespace RAPTOR