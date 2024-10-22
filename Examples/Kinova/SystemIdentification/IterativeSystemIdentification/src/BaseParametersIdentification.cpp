#include "BaseParametersIdentification.h"

namespace RAPTOR {

// // constructor
// BaseParametersIdentification::BaseParametersIdentification()
// {
// }

// // destructors
// BaseParametersIdentification::~BaseParametersIdentification()
// {
// }

bool BaseParametersIdentification::set_parameters(
    const Model& model_input,
    const std::shared_ptr<MatX>& posPtr_input,
    const std::shared_ptr<MatX>& velPtr_input,
    const std::shared_ptr<MatX>& accPtr_input,
    const std::shared_ptr<MatX>& torquePtr_input,
    std::shared_ptr<QRDecompositionSolver> regroupPtr_input,
    const bool include_offset_input
)
{ 
    enable_hessian = false;

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

    regroupPtr_ = regroupPtr_input;
    include_offset = include_offset_input;

    // initialize observation matrices
    FullObservationMatrix = MatX::Zero(N * Nact, 10 * Nact);
    for (int i = 0; i < N; i++) {
        const VecX& q = posPtr_->col(i);
        const VecX& v = velPtr_->col(i);
        const VecX& a = accPtr_->col(i);

        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, 
            q, v, a);

        FullObservationMatrix.middleRows(i * Nact, Nact) = 
            dataPtr_->jointTorqueRegressor;
    }

    // Perform regrouping (assume regroupPtr_ has been initialized)
    if (regroupPtr_->RegroupMatrix.rows() != FullObservationMatrix.cols()) {
        throw std::invalid_argument("Error: QRDecompositionSolver not initialized properly!");
    }

    RegroupedObservationMatrix = 
        FullObservationMatrix * regroupPtr_->RegroupMatrix;

    // directly set up initial condition here
    int n = regroupPtr_->dim_id + // base parameters 
            regroupPtr_->dim_d +  // dependent parameters
            3 * Nact;             // friction, damping, armature
    if (include_offset) {
        n += Nact; // offset
    }

    x0 = VecX::Zero(n);

        // initial guess for independent parameters
        // is just what is included in the original urdf
    x0.head(regroupPtr_->dim_id) = regroupPtr_->beta;
    x0.segment(regroupPtr_->dim_id, regroupPtr_->dim_d) = regroupPtr_->phi_d;

        // initial guess for motor friction parameters is just 0

    // initialize LMI constraints for all links
    constraintsPtrVec_.push_back(std::make_unique<RegroupedLMIConstraints>(regroupPtr_input,
                                                                           modelPtr_->nv,
                                                                           n));       
    constraintsNameVec_.push_back("Regrouped LMI constraints"); 

    return true;
}

bool BaseParametersIdentification::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n = regroupPtr_->dim_id + // base parameters 
        regroupPtr_->dim_d +  // dependent parameters
        3 * Nact;             // friction, damping, armature
    if (include_offset) {
        n += Nact;
    }
    numVars= n;

    // number of inequality constraint
    numCons = 0;
    for ( Index i = 0; i < constraintsPtrVec_.size(); i++ ) {
        numCons += constraintsPtrVec_[i]->m;
    }
    m = numCons;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool BaseParametersIdentification::get_bounds_info(
    Index n, 
    Number* x_l, 
    Number* x_u,
    Index m, 
    Number* g_l, 
    Number* g_u
)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_bounds_info!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in get_bounds_info!");
    }

    // use base class function to set bounds in g_l and g_u
    Optimizer::get_bounds_info(n, x_l, x_u, m, g_l, g_u);

    // set variable bounds (overwrite previous bounds in x_l and x_u)
        // independent inertial parameters (after regrouping)
    for (Index i = 0; i < regroupPtr_->dim_id; i++) {
        if (regroupPtr_->beta(i) > 0) {
            x_l[i] = (1 - default_maximum_uncertainty) * regroupPtr_->beta(i);
            x_u[i] = (1 + default_maximum_uncertainty) * regroupPtr_->beta(i);
        }
        else {
            x_l[i] = (1 + default_maximum_uncertainty) * regroupPtr_->beta(i);
            x_u[i] = (1 - default_maximum_uncertainty) * regroupPtr_->beta(i);
        }
    }

        // dependent inertial parameters (after regrouping)
    for (Index i = 0; i < regroupPtr_->dim_d; i++) {
        if (regroupPtr_->phi_d(i) > 0) {
            x_l[i] = (1 - default_maximum_uncertainty) * regroupPtr_->phi_d(i);
            x_u[i] = (1 + default_maximum_uncertainty) * regroupPtr_->phi_d(i);
        }
        else {
            x_l[i] = (1 + default_maximum_uncertainty) * regroupPtr_->phi_d(i);
            x_u[i] = (1 - default_maximum_uncertainty) * regroupPtr_->phi_d(i);
        }
    }

        // static friction >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[regroupPtr_->dim_id + i] = 0.0;
        x_u[regroupPtr_->dim_id + i] = 50.0;
    }

        // damping >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[regroupPtr_->dim_id + Nact + i] = 0.0;
        x_u[regroupPtr_->dim_id + Nact + i] = 50.0;
    }

        // armature >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[regroupPtr_->dim_id + 2 * Nact + i] = 0.0;
        x_u[regroupPtr_->dim_id + 2 * Nact + i] = 50.0;
    }

        // static friction offset
    if (include_offset) {
        for (Index i = 0; i < Nact; i++) {
            x_l[regroupPtr_->dim_id + 3 * Nact + i] = -50.0;
            x_u[regroupPtr_->dim_id + 3 * Nact + i] = 50.0;
        }
    }

    return true;
}

bool BaseParametersIdentification::eval_f(
    Index n,
    const Number *x,
    bool new_x,
    Number &obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    const int& dim_id = regroupPtr_->dim_id;
    const int& dim_d = regroupPtr_->dim_d;

    const VecX& beta = z.head(dim_id);
    // const VecX& phi_d = z.segment(dim_id, Nact);
    const VecX& friction = z.segment(dim_id + dim_d, Nact);
    const VecX& damping = z.segment(dim_id + dim_d + Nact, Nact);
    const VecX& armature = z.segment(dim_id + dim_d + 2 * Nact, Nact);
    VecX offset = VecX::Zero(Nact);
    if (include_offset) {
        offset = z.tail(Nact);
    }

    VecX tau_inertials = RegroupedObservationMatrix * beta;

    obj_value = 0;

    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);
        const VecX& tau = torquePtr_->col(i);
        const VecX& tau_inertial = tau_inertials.segment(i * Nact, Nact);

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

bool BaseParametersIdentification::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);

    const int& dim_id = regroupPtr_->dim_id;
    const int& dim_d = regroupPtr_->dim_d;

    const VecX& beta = z.head(dim_id);
    // const VecX& phi_d = z.segment(dim_id, Nact);
    const VecX& friction = z.segment(dim_id + dim_d, Nact);
    const VecX& damping = z.segment(dim_id + dim_d + Nact, Nact);
    const VecX& armature = z.segment(dim_id + dim_d + 2 * Nact, Nact);
    VecX offset = VecX::Zero(Nact);
    if (include_offset) {
        offset = z.tail(Nact);
    }

    VecX tau_inertials = RegroupedObservationMatrix * beta;

    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);
        const VecX& tau = torquePtr_->col(i);
        const VecX& tau_inertial = tau_inertials.segment(i * Nact, Nact);

        VecX total_friction_force = 
            friction.cwiseProduct(q_d.cwiseSign()) +
            damping.cwiseProduct(q_d) +
            armature.cwiseProduct(q_dd) +
            offset;

        VecX tau_estimated = tau_inertial + total_friction_force;
        VecX tau_diff = tau_estimated - tau;

        grad_f_vec.head(dim_id) += tau_diff.transpose() * (RegroupedObservationMatrix.middleRows(i * Nact, Nact));
        grad_f_vec.segment(dim_id + dim_d, Nact) += tau_diff.cwiseProduct(q_d.cwiseSign());
        grad_f_vec.segment(dim_id + dim_d + Nact, Nact) += tau_diff.cwiseProduct(q_d);
        grad_f_vec.segment(dim_id + dim_d + 2 * Nact, Nact) += tau_diff.cwiseProduct(q_dd);
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
bool BaseParametersIdentification::eval_hess_f(
   Index         n,
   const Number* x,
   bool          new_x,
   MatX&         hess_f
)
{
    if (n != numVars) {
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    const int& dim_id = regroupPtr_->dim_id;
    const int& dim_d = regroupPtr_->dim_d;
    const int dim = dim_id + dim_d;

    hess_f = MatX::Zero(n, n);

    for (Index i = 0; i < N; i++) {
        const VecX& q_d = velPtr_->col(i);
        const VecX& q_dd = accPtr_->col(i);
        const MatX& romi = RegroupedObservationMatrix.middleRows(i * Nact, Nact);
        const MatX romi_T = romi.transpose();
        
        hess_f.topLeftCorner(dim_id, dim_id) += 
            romi_T * romi;
        hess_f.block(0, dim, dim_id, Nact).array() += 
            (romi.array().colwise() * q_d.cwiseSign().array()).transpose();
        hess_f.block(0, dim + Nact, dim_id, Nact).array() +=
            (romi.array().colwise() * q_d.array()).transpose();
        hess_f.block(0, dim + 2 * Nact, dim_id, Nact).array() +=
            (romi.array().colwise() * q_dd.array()).transpose();
        if (include_offset) {
            hess_f.block(0, dim + 3 * Nact, dim_id, Nact).array() +=
                romi.array().transpose();
        }

        hess_f.diagonal().segment(dim, Nact) += q_d.cwiseSign().cwiseAbs2();
        hess_f.block(dim, dim + Nact, Nact, Nact).diagonal().array() += q_d.cwiseSign().array() * q_d.array();
        hess_f.block(dim, dim + 2 * Nact, Nact, Nact).diagonal().array() += q_d.cwiseSign().array() * q_dd.array();
        if (include_offset) {
            hess_f.block(dim, dim + 3 * Nact, Nact, Nact).diagonal() += q_d.cwiseSign();
        }

        hess_f.diagonal().segment(dim + Nact, Nact) += q_d.cwiseAbs2();
        hess_f.block(dim + Nact, dim + 2 * Nact, Nact, Nact).diagonal().array() += q_d.array() * q_dd.array();
        if (include_offset) {
            hess_f.block(dim + Nact, dim + 3 * Nact, Nact, Nact).diagonal() += q_d;
        }

        hess_f.diagonal().segment(dim + 2 * Nact, Nact) += q_dd.cwiseAbs2();
        if (include_offset) {
            hess_f.block(dim + 2 * Nact, dim + 3 * Nact, Nact, Nact).diagonal() += q_dd;
        }

        if (include_offset) {
            hess_f.diagonal().tail(Nact) += VecX::Ones(Nact);
        }
    }

    return true;
}
// [TNLP_eval_hess_f]

}; // namespace RAPTOR
