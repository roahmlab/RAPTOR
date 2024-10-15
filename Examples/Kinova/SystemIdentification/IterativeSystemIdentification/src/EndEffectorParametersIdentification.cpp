#include "EndEffectorParametersIdentification.h"

namespace RAPTOR {

// // constructor
// EndEffectorParametersIdentification::EndEffectorParametersIdentification()
// {
// }

// // destructors
// EndEffectorParametersIdentification::~EndEffectorParametersIdentification()
// {
// }

bool EndEffectorParametersIdentification::set_parameters(
    const Model& model_input,
    const std::shared_ptr<MatX>& posPtr_input,
    const std::shared_ptr<MatX>& velPtr_input,
    const std::shared_ptr<MatX>& accPtr_input,
    const std::shared_ptr<MatX>& torquePtr_input,
    const std::shared_ptr<VecX>& phiPtr_input,
    const bool include_offset_input
)
{ 
    enable_hessian = false;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    // only choose optimize the end-effector
    Nact = modelPtr_->nv;

    phiPtr_ = phiPtr_input;

    if (phiPtr_->size() < 12 * Nact) {
        throw std::invalid_argument("Error: phiPtr_ size is too small.");
    }

    friction = phiPtr_->segment(10 * Nact , Nact);
    damping = phiPtr_->segment(11 * Nact, Nact);
    armature = phiPtr_->segment(12 * Nact, Nact);
    offset = VecX::Zero(Nact);
    if (include_offset) {
        offset = phiPtr_->tail(Nact);
    }


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

    // directly set up initial condition here
    int n = 10;  // End-effector parameters

    x0 = VecX::Zero(n);

    // initial guess for independent parameters
    // is just what is included in the optimised in the full parameters optimise
    x0 = phiPtr_->segment(10 * (Nact-1), 10);
    std::cout << "end_effector parameters"<< x0.transpose() <<std::endl;

    // initial guess for motor friction parameters is just 0

    // initialize LMI constraints for all links
    constraintsPtrVec_.push_back(std::make_unique<LMIConstraints>(1,
                                                                  n));       
    constraintsNameVec_.push_back("LMI constraints"); 

    return true;
}

bool EndEffectorParametersIdentification::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n =10;       // End-effector parameters 
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

bool EndEffectorParametersIdentification::get_bounds_info(
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
    // inertial parameters
    for (Index i = 0; i < 10; i++) {
        if (x0(i) > 0) {
            x_l[i] = (1 - default_maximum_uncertainty) * x0(i);
            x_u[i] = (1 + default_maximum_uncertainty) * x0(i);
        }
        else {
            x_l[i] = (1 + default_maximum_uncertainty) * x0(i);
            x_u[i] = (1 - default_maximum_uncertainty) * x0(i);
        }
    }

    return true;
}

bool EndEffectorParametersIdentification::eval_f(
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

    phiPtr_->segment((Nact -1)* 10, 10)= z;

    // the tau_inertials of end+effector
    VecX tau_inertials = FullObservationMatrix * phiPtr_->head(10 * Nact);


    obj_value = 0;

    for (Index i =0 ; i < N; i++) {
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

bool EndEffectorParametersIdentification::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);

    phiPtr_->segment(10 * (Nact-1), 10)= z;

    VecX tau_inertials =FullObservationMatrix * phiPtr_->head(10 * Nact);

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

        grad_f_vec.head(10) += tau_diff.transpose() * (FullObservationMatrix.block(i * Nact, 
                                                                                   10* (Nact -1),
                                                                                   Nact,
                                                                                   10));

    }

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i);
    }

    return true;
}

// [TNLP_eval_hess_f]
// return the hessian of the objective function hess_{x} f(x) as a dense matrix
bool EndEffectorParametersIdentification::eval_hess_f(
   Index         n,
   const Number* x,
   bool          new_x,
   MatX&         hess_f
)
{
    if (n != numVars) {
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    // const int& dim_id = regroupPtr_->dim_id;
    // const int& dim_d = regroupPtr_->dim_d;
    // const int dim = dim_id + dim_d;

    hess_f = MatX::Zero(n, n);  

    // for (Index i = 0; i < N; i++) {
    //     const VecX& q_d = velPtr_->col(i);
    //     const VecX& q_dd = accPtr_->col(i);
    //     const MatX& romi = FullObservationMatrix.block( i * Nact, 
    //                                                     10* (Nact -1),
    //                                                     Nact,
    //                                                     10);
    //     const MatX romi_T = romi.transpose();
        
    //     hess_f += 
    //         romi_T * romi;
  
    // }

    const MatX& romi = FullObservationMatrix.block( 0, 
                                                        10* (Nact -1),
                                                        FullObservationMatrix.rows(),
                                                        10);
    const MatX romi_T = romi.transpose();
    
    hess_f = 
        romi_T * romi;

    return true;
}
// [TNLP_eval_hess_f]

}; // namespace RAPTOR
