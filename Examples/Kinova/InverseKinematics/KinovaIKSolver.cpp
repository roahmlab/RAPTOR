#include "KinovaIKSolver.h"

namespace RAPTOR {
namespace Kinova {

// // constructor
// KinovaIKSolver::KinovaIKSolver()
// {
// }


// // destructor
// KinovaIKSolver::~KinovaIKSolver()
// {
// }

// [TNLP_set_parameters]
bool KinovaIKSolver::set_parameters(
    const VecX& x0_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    const Transform& desiredTransform_input
)
{
    enable_hessian = true;
    x0 = x0_input;

    trajPtr_ = std::make_shared<Plain>(model_input.nq);
    
    // read joint limits from KinovaConstants.h
    VecX JOINT_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS);
    VecX JOINT_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS);

    // Joint limits
    constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
                                                               JOINT_LIMITS_LOWER_VEC, 
                                                               JOINT_LIMITS_UPPER_VEC));
    constraintsNameVec_.push_back("joint limits");

    // End effector kinematics constraints
    constraintsPtrVec_.push_back(std::make_unique<KinematicsConstraints>(trajPtr_,
                                                                         &model_input,
                                                                         jtype_input,
                                                                         model_input.nq, // the last joint
                                                                         0,
                                                                         desiredTransform_input));   
    constraintsNameVec_.push_back("kinematics constraints");                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                                                                        
    assert(x0.size() == trajPtr_->varLength);

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
// returns some info about the nlp
bool KinovaIKSolver::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // number of decision variables
    numVars = trajPtr_->varLength;
    n = numVars;

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
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool KinovaIKSolver::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
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

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = JOINT_LIMITS_LOWER[i];
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = JOINT_LIMITS_UPPER[i];
    }

    if (constraintsPtrVec_.size() != constraintsNameVec_.size()) {
        THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ and constraintsNameVec_ have different sizes!");
    }

    // compute bounds for all constraints
    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        try {
            constraintsPtrVec_[c]->compute_bounds();
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in get_bounds_info! Check previous error message.");
        }

        if (constraintsPtrVec_[c]->m != constraintsPtrVec_[c]->g_lb.size() || 
            constraintsPtrVec_[c]->m != constraintsPtrVec_[c]->g_ub.size()) {
            THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ have different sizes!");
        }

        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g_l[iter] = constraintsPtrVec_[c]->g_lb(i);
            g_u[iter] = constraintsPtrVec_[c]->g_ub(i);
            iter++;
        }
    }

    // report constraints distribution
    std::cout << "Dimension of each constraints and their locations: \n";
    iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        std::cout << constraintsNameVec_[c] << ": ";
        std::cout << constraintsPtrVec_[c]->m << " [";
        std::cout << iter << " ";
        iter += constraintsPtrVec_[c]->m;
        std::cout << iter << "]\n";
    }

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_eval_f]
// returns the value of the objective function
bool KinovaIKSolver::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, false);

    const VecX& q = trajPtr_->q(0);

    obj_value = 0.5 * q.dot(q);

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool KinovaIKSolver::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, true);

    const VecX& q = trajPtr_->q(0);
    const MatX& pq_pz = trajPtr_->pq_pz(0);

    VecX grad = q.transpose() * pq_pz;
    for(Index i = 0; i < n; i++){
        grad_f[i] = grad(i);
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_hess_f]
// return the hessian of the objective function hess_{x} f(x) as a dense matrix
bool KinovaIKSolver::eval_hess_f(
   Index         n,
   const Number* x,
   bool          new_x,
   MatX&         hess_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, true, true);

    const VecX& q = trajPtr_->q(0);
    const MatX& pq_pz = trajPtr_->pq_pz(0);
    // we know the following is just 0
    // const Eigen::Array<MatX, Eigen::Dynamic, 1>& pq_pz_pz = trajPtr_->pq_pz_pz.col(0);

    hess_f = pq_pz.transpose() * pq_pz;

    return true;
}
// [TNLP_eval_hess_f]

}; // namespace Kinova
}; // namespace RAPTOR