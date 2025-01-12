#include "KinovaLongerHorizonOptimizer.h"

namespace RAPTOR {
namespace Kinova {

// // constructor
// KinovaLongerHorizonOptimizer::KinovaLongerHorizonOptimizer()
// {
// }


// // destructor
// KinovaLongerHorizonOptimizer::~KinovaLongerHorizonOptimizer()
// {
// }

// [TNLP_set_parameters]
bool KinovaLongerHorizonOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const int degree_input,
    const Model& model_input, 
    const VecX& q0_input,
    const VecX& qT_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input,
    const VecX& joint_limits_buffer_input,
    const VecX& velocity_limits_buffer_input,
    const VecX& torque_limits_buffer_input,
    const bool include_gripper_or_not,
    const double collision_buffer_input
 ) 
{
    enable_hessian = true;
    x0 = x0_input;

    trajPtr_ = std::make_shared<PiecewiseBezierCurves>(T_input, 
                                                       N_input, 
                                                       model_input.nq, 
                                                       degree_input, 
                                                       q0_input,
                                                       qT_input);

    idPtr_ = std::make_shared<InverseDynamics>(model_input,
                                               trajPtr_);
    
    // read joint limits from KinovaConstants.h
    VecX JOINT_LIMITS_LOWER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS)) + 
        joint_limits_buffer_input;

    VecX JOINT_LIMITS_UPPER_VEC =
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS)) -
        joint_limits_buffer_input;

    // read velocity limits from KinovaConstants.h
    VecX VELOCITY_LIMITS_LOWER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(VELOCITY_LIMITS_LOWER, NUM_JOINTS)) +
        velocity_limits_buffer_input;

    VecX VELOCITY_LIMITS_UPPER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(VELOCITY_LIMITS_UPPER, NUM_JOINTS)) -
        velocity_limits_buffer_input;

    // read torque limits from KinovaConstants.h
    VecX TORQUE_LIMITS_LOWER_VEC = 
        Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_LOWER, NUM_JOINTS) + 
        torque_limits_buffer_input;

    VecX TORQUE_LIMITS_UPPER_VEC = 
        Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_UPPER, NUM_JOINTS) -
        torque_limits_buffer_input;

    // // Joint limits
    // constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
    //                                                            JOINT_LIMITS_LOWER_VEC, 
    //                                                            JOINT_LIMITS_UPPER_VEC));
    // constraintsNameVec_.push_back("joint limits");

    // // Velocity limits
    // constraintsPtrVec_.push_back(std::make_unique<VelocityLimits>(trajPtr_, 
    //                                                               VELOCITY_LIMITS_LOWER_VEC, 
    //                                                               VELOCITY_LIMITS_UPPER_VEC));
    // constraintsNameVec_.push_back("velocity limits");        

    // // Torque limits
    // constraintsPtrVec_.push_back(std::make_unique<TorqueLimits>(trajPtr_, 
    //                                                             idPtr_,
    //                                                             TORQUE_LIMITS_LOWER_VEC, 
    //                                                             TORQUE_LIMITS_UPPER_VEC));
    // constraintsNameVec_.push_back("torque limits");                                                            

    // Customized constraints (collision avoidance with obstacles)
    constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_,
                                                                               model_input,
                                                                               boxCenters_input,
                                                                               boxOrientation_input,
                                                                               boxSize_input,
                                                                               include_gripper_or_not,
                                                                               collision_buffer_input));   
    constraintsNameVec_.push_back("obstacle avoidance constraints");

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
// returns some info about the nlp
bool KinovaLongerHorizonOptimizer::get_nlp_info(
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

// [TNLP_eval_f]
// returns the value of the objective function
bool KinovaLongerHorizonOptimizer::eval_f(
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

    obj_value = 0.0;

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool KinovaLongerHorizonOptimizer::eval_grad_f(
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
    
    for(Index i = 0; i < n; i++){
        grad_f[i] = 0.0;
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_hess_f]
// return the hessian of the objective function hess_{x} f(x) as a dense matrix
bool KinovaLongerHorizonOptimizer::eval_hess_f(
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

    hess_f = MatX::Zero(n, n);

    return true;
}
// [TNLP_eval_hess_f]

}; // namespace Kinova
}; // namespace RAPTOR