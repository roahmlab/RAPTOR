#include "DigitSystemIdentification.h"

namespace RAPTOR {
// namespace DigitWholeBodySysID {
namespace Digit {

// // constructor
// DigitSystemIdentification::DigitSystemIdentification()
// {
// }

// // destructors
// DigitSystemIdentification::~DigitSystemIdentification()
// {
// }

bool DigitSystemIdentification::set_parameters(
    const Model& model_input,
    const std::shared_ptr<MatX>& posPtr_input,
    const std::shared_ptr<MatX>& velPtr_input,
    const std::shared_ptr<MatX>& accPtr_input,
    const std::shared_ptr<MatX>& torquePtr_input,
    const char stanceLeg_input,
    const Transform& stance_foot_T_des_input
)
{ 
    enable_hessian = false;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    // ddcPtr_ = std::make_shared<DigitWholeBodyDynamicsConstraints>(modelPtr_, 
    //                                                               stanceLeg_input,
    //                                                               stance_foot_T_des_input);
    ddcPtr_ = std::make_shared<DigitDynamicsConstraints>(modelPtr_, 
                                                         stanceLeg_input,
                                                         stance_foot_T_des_input);

    Nact = NUM_INDEPENDENT_JOINTS;

    posPtr_ = posPtr_input;
    velPtr_ = velPtr_input;
    accPtr_ = accPtr_input;
    torquePtr_ = torquePtr_input;

    if (posPtr_->rows() != velPtr_->rows() || 
        posPtr_->rows() != accPtr_->rows() || 
        posPtr_->rows() != modelPtr_->nv) {
        throw std::invalid_argument("Error: input data matrices have different number of rows");
    }

    if (torquePtr_->rows() != Nact) {
        throw std::invalid_argument("Error: torque data matrix has different number of rows");
    }

    if (posPtr_->cols() != velPtr_->cols() || 
        posPtr_->cols() != accPtr_->cols() || 
        posPtr_->cols() != torquePtr_->cols()) {
        throw std::invalid_argument("Error: input data matrices have different number of columns");
    }

    N = posPtr_->cols();

    // find links that are not trivial (non-zero mass and inertia)
    nontrivialLinkIds.clear();
    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1;
        if (modelPtr_->inertias[pinocchio_joint_id].mass() > 0.0) {
            nontrivialLinkIds.push_back(i);
        }
    }

    // directly set up initial condition here
    int n = 10 * nontrivialLinkIds.size() + 3 * Nact;

    x0 = VecX::Zero(n);
    for (int i = 0; i < nontrivialLinkIds.size(); i++) {
        const int pinocchio_joint_id = nontrivialLinkIds[i] + 1;
        x0.segment(10 * i, 10) = modelPtr_->inertias[pinocchio_joint_id].toDynamicParameters();
    }
    x0.segment(10 * nontrivialLinkIds.size(), Nact) = ddcPtr_->get_independent_vector(modelPtr_->friction);
    x0.segment(10 * nontrivialLinkIds.size() + Nact, Nact) = ddcPtr_->get_independent_vector(modelPtr_->damping);
    x0.segment(10 * nontrivialLinkIds.size() + 2 * Nact, Nact) = ddcPtr_->get_independent_vector(modelPtr_->armature);

    // initial guess for motor friction parameters is just 0

    // initialize observation matrices
    FullObservationMatrix = MatX::Zero(N * Nact, 13 * modelPtr_->nv);

    for (int i = 0; i < N; i++) {
        MatX SubObservationMatrix = MatX::Zero(modelPtr_->nv, 13 * modelPtr_->nv);

        VecX q = posPtr_->col(i);
        VecX v = velPtr_->col(i);
        VecX a = accPtr_->col(i);

        ddcPtr_->setupJointPositionVelocityAcceleration(q, v, a, false);

        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, 
            q, v, a);

        SubObservationMatrix.leftCols(10 * modelPtr_->nv) = dataPtr_->jointTorqueRegressor;
        SubObservationMatrix.middleCols(10 * modelPtr_->nv, modelPtr_->nv).diagonal() = v.cwiseSign();
        SubObservationMatrix.middleCols(11 * modelPtr_->nv, modelPtr_->nv).diagonal() = v;
        SubObservationMatrix.middleCols(12 * modelPtr_->nv, modelPtr_->nv).diagonal() = a;

        MatX IndependentObservationMatrix(NUM_INDEPENDENT_JOINTS, 13 * modelPtr_->nv);
        MatX DependentObservationMatrix(NUM_DEPENDENT_JOINTS, 13 * modelPtr_->nv);
        ddcPtr_->get_independent_rows(IndependentObservationMatrix, SubObservationMatrix);
        ddcPtr_->get_dependent_rows(DependentObservationMatrix, SubObservationMatrix);
        
        FullObservationMatrix.middleRows(i * Nact, Nact) = 
            IndependentObservationMatrix +
            ddcPtr_->P_dep.transpose() * DependentObservationMatrix;
    }

    // initialize LMI constraints for all links
    constraintsPtrVec_.push_back(std::make_unique<LMIConstraints>(nontrivialLinkIds.size(), n));       
    constraintsNameVec_.push_back("LMI constraints"); 

    return true;
}

bool DigitSystemIdentification::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n = 10 * nontrivialLinkIds.size() + 3 * Nact;
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

bool DigitSystemIdentification::get_bounds_info(
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
    for (Index i = 0; i < 10 * nontrivialLinkIds.size(); i++) {
        if (x0(i) > 0) {
            x_l[i] = (1 - default_maximum_uncertainty) * x0(i);
            x_u[i] = (1 + default_maximum_uncertainty) * x0(i);
        }
        else {
            x_l[i] = (1 + default_maximum_uncertainty) * x0(i);
            x_u[i] = (1 - default_maximum_uncertainty) * x0(i);
        }
    }

        // static friction >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[10 * nontrivialLinkIds.size() + i] = 0;
        x_u[10 * nontrivialLinkIds.size() + i] = 50.0;
    }

        // damping >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[10 * nontrivialLinkIds.size() + Nact + i] = 0;
        x_u[10 * nontrivialLinkIds.size() + Nact + i] = 50.0;
    }

        // armature >= 0
    for (Index i = 0; i < Nact; i++) {
        x_l[10 * nontrivialLinkIds.size() + 2 * Nact + i] = 0;
        x_u[10 * nontrivialLinkIds.size() + 2 * Nact + i] = 50.0;
    }

    return true;
}

bool DigitSystemIdentification::eval_f(
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

    VecX phi_full = VecX::Zero(13 * modelPtr_->nv);
    for (Index i = 0; i < nontrivialLinkIds.size(); i++) {
        phi_full.segment(10 * nontrivialLinkIds[i], 10) = z.segment(10 * i, 10);
    }
    VecX friction_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(friction_full, 
                                     z.segment(10 * nontrivialLinkIds.size(), Nact));
    phi_full.segment(10 * modelPtr_->nv, modelPtr_->nv) = friction_full;
    VecX damping_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(damping_full, 
                                     z.segment(10 * nontrivialLinkIds.size() + Nact, Nact));
    phi_full.segment(11 * modelPtr_->nv, modelPtr_->nv) = damping_full;
    VecX armature_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(armature_full, 
                                     z.segment(10 * nontrivialLinkIds.size() + 2 * Nact, Nact));
    phi_full.segment(12 * modelPtr_->nv, modelPtr_->nv) = armature_full;

    const VecX tau_estimated = FullObservationMatrix * phi_full;

    obj_value = 0;

    for (Index i = 0; i < N; i++) {
        const VecX& tau = torquePtr_->col(i);
        obj_value += 0.5 * (tau_estimated.segment(i * Nact, Nact) - tau).squaredNorm();
    } 

    obj_value /= N;

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool DigitSystemIdentification::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);

    VecX phi_full = VecX::Zero(13 * modelPtr_->nv);
    for (Index i = 0; i < nontrivialLinkIds.size(); i++) {
        phi_full.segment(10 * nontrivialLinkIds[i], 10) = z.segment(10 * i, 10);
    }
    VecX friction_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(friction_full, 
                                     z.segment(10 * nontrivialLinkIds.size(), Nact));
    phi_full.segment(10 * modelPtr_->nv, modelPtr_->nv) = friction_full;
    VecX damping_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(damping_full, 
                                     z.segment(10 * nontrivialLinkIds.size() + Nact, Nact));
    phi_full.segment(11 * modelPtr_->nv, modelPtr_->nv) = damping_full;
    VecX armature_full(modelPtr_->nv);
    ddcPtr_->fill_independent_vector(armature_full, 
                                     z.segment(10 * nontrivialLinkIds.size() + 2 * Nact, Nact));
    phi_full.segment(12 * modelPtr_->nv, modelPtr_->nv) = armature_full;

    const VecX tau_estimated = FullObservationMatrix * phi_full;

    for (Index i = 0; i < N; i++) {
        const VecX& tau = torquePtr_->col(i);

        const VecX grad_f_vec_full = 
            (tau_estimated.segment(i * Nact, Nact) - tau).transpose() * 
                      FullObservationMatrix.middleRows(i * Nact, Nact);
        
        for (Index j = 0; j < nontrivialLinkIds.size(); j++) {
            grad_f_vec.segment(10 * j, 10) += grad_f_vec_full.segment(10 * nontrivialLinkIds[j], 10);
        }
        grad_f_vec.segment(10 * nontrivialLinkIds.size(), Nact) +=
            ddcPtr_->get_independent_vector(grad_f_vec_full.segment(10 * modelPtr_->nv, modelPtr_->nv));
        grad_f_vec.segment(10 * nontrivialLinkIds.size() + Nact, Nact) +=
            ddcPtr_->get_independent_vector(grad_f_vec_full.segment(11 * modelPtr_->nv, modelPtr_->nv));
        grad_f_vec.segment(10 * nontrivialLinkIds.size() + 2 * Nact, Nact) +=
            ddcPtr_->get_independent_vector(grad_f_vec_full.segment(12 * modelPtr_->nv, modelPtr_->nv));
    }

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i) / N;
    }

    return true;
}

}; // namespace DigitWholeBodySysID
}; // namespace RAPTOR
