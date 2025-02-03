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
    const Transform& desiredTransform_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input,
    const Transform endT_input,
    const bool include_gripper_or_not,
    const double collision_buffer_input,
    Eigen::VectorXi jtype_input
)
{
    enable_hessian = true;
    x0 = x0_input;

    trajPtr_ = std::make_shared<Plain>(model_input.nq);
    
    // read joint limits from KinovaConstants.h
    VecX JOINT_LIMITS_LOWER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS));

    VecX JOINT_LIMITS_UPPER_VEC =
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS));

    // In this instance, the decision variables are the joint angles of the robot directly.
    // You can see that we inherited the get_bounds_info function from the Optimizer class
    // and we are modifing it (x_lb and x_ub) to directly set the bounds for the joint limits.
    // As a result, the following joint limits constraints are not needed.
    // It is equivalent to turn on this constraint and remove the inherited get_bounds_info function.

    // // Joint limits
    // constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
    //                                                            JOINT_LIMITS_LOWER_VEC, 
    //                                                            JOINT_LIMITS_UPPER_VEC));
    // constraintsNameVec_.push_back("joint limits");

    // End effector kinematics constraints
    constraintsPtrVec_.push_back(std::make_unique<KinematicsConstraints>(trajPtr_,
                                                                         &model_input,
                                                                         model_input.nq, // the last joint
                                                                         0,
                                                                         desiredTransform_input,
                                                                         endT_input,
                                                                         jtype_input));   
    constraintsNameVec_.push_back("kinematics constraints");   

    if (boxCenters_input.size() != boxOrientation_input.size() || 
        boxCenters_input.size() != boxSize_input.size()) {
        throw std::invalid_argument("boxCenters_input, boxOrientation_input, and boxSize_input have different sizes!");
    }                      

    if (boxCenters_input.size() > 0) {
        // Collision avoidance constraints
        constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_,
                                                                                   model_input,
                                                                                   boxCenters_input,
                                                                                   boxOrientation_input,
                                                                                   boxSize_input,
                                                                                   include_gripper_or_not,
                                                                                   collision_buffer_input,
                                                                                   jtype_input));
        constraintsNameVec_.push_back("collision avoidance constraints");
    }
                                                                                                                                                                                                                                                                                                                                                                        
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

    // number of constraints
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

    // In this instance, the decision variables are the joint angles of the robot directly.
    // You can see that we inherited the get_bounds_info function from the Optimizer class
    // and we are modifing it (x_lb and x_ub) to directly set the bounds for the joint limits.
    // As a result, the following joint limits constraints are not needed.
    // It is equivalent to turn on this constraint and remove the inherited get_bounds_info function.

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
    if (display_info) {
        std::cout << "Dimension of each constraints and their locations: \n";
        iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            std::cout << constraintsNameVec_[c] << ": ";
            std::cout << constraintsPtrVec_[c]->m << " [";
            std::cout << iter << " ";
            iter += constraintsPtrVec_[c]->m;
            std::cout << iter << "]\n";
        }
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

    // minimize the distance to a given initial configuration
    // kinova has 4 continuous joints
    obj_value = pow(Utils::wrapToPi(q(0) - x0(0)), 2) + // These are continuous joints
                pow(Utils::wrapToPi(q(2) - x0(2)), 2) + 
                pow(Utils::wrapToPi(q(4) - x0(4)), 2) + 
                pow(Utils::wrapToPi(q(6) - x0(6)), 2) + 
                pow(q(1) - x0(1), 2) +                  // These are not continuous joints
                pow(q(3) - x0(3), 2) + 
                pow(q(5) - x0(5), 2);
    obj_value = 0.5 * obj_value;

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

    for (Index i = 0; i < n; i++) {
        // kinova has 4 continuous joints
        if (i % 2 == 0) {
            grad_f[i] = Utils::wrapToPi(q(i) - x0(i)) * pq_pz(i, i);
        }
        else {
            grad_f[i] = (q(i) - x0(i)) * pq_pz(i, i);
        }
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