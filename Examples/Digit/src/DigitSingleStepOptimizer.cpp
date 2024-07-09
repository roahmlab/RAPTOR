#include "DigitSingleStepOptimizer.h"

namespace IDTO {
namespace Digit {

// // constructor
// DigitSingleStepOptimizer::DigitSingleStepOptimizer()
// {
// }


// // destructor
// DigitSingleStepOptimizer::~DigitSingleStepOptimizer()
// {
// }

bool DigitSingleStepOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const TimeDiscretization time_discretization_input,
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    const GaitParameters& gp_input
 ) 
{
    x0 = x0_input;
                                          
    trajPtr_ = std::make_shared<BezierCurves>(T_input, 
                                              N_input, 
                                              NUM_INDEPENDENT_JOINTS, 
                                              time_discretization_input, 
                                              degree_input);                                   
    
    // add v_reset and lambda_reset to the end of the decision variables                                         
    trajPtr_->varLength += NUM_JOINTS + NUM_DEPENDENT_JOINTS;
    
    // stance foot is left foot by default
    char stanceLeg = 'L';
    Transform stance_foot_T_des(3, -M_PI / 2);
    dcidPtr_ = std::make_shared<DigitConstrainedInverseDynamics>(model_input, 
                                                                 trajPtr_,
                                                                 NUM_DEPENDENT_JOINTS, 
                                                                 jtype_input, 
                                                                 stanceLeg, 
                                                                 stance_foot_T_des);                                                          
    cidPtr_ = dcidPtr_; // convert to base class
    
    // convert joint limits from degree to radian
    VecX JOINT_LIMITS_LOWER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_LOWER_VEC(i) = Utils::deg2rad(JOINT_LIMITS_LOWER[i]);
    }

    // convert joint limits from degree to radian   
    VecX JOINT_LIMITS_UPPER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_UPPER_VEC(i) = Utils::deg2rad(JOINT_LIMITS_UPPER[i]);
    }

    VecX TORQUE_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_LOWER, NUM_INDEPENDENT_JOINTS);
    VecX TORQUE_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_UPPER, NUM_INDEPENDENT_JOINTS);

    constraintsPtrVec_.clear();  
     
    // Torque limits
    constraintsPtrVec_.push_back(std::make_unique<TorqueLimits>(trajPtr_, 
                                                                cidPtr_, 
                                                                TORQUE_LIMITS_LOWER_VEC, 
                                                                TORQUE_LIMITS_UPPER_VEC));        
    constraintsNameVec_.push_back("torque limits");

    // Joint limits
    constraintsPtrVec_.push_back(std::make_unique<ConstrainedJointLimits>(trajPtr_, 
                                                                          cidPtr_->dcPtr_, 
                                                                          JOINT_LIMITS_LOWER_VEC, 
                                                                          JOINT_LIMITS_UPPER_VEC));      
    constraintsNameVec_.push_back("joint limits");                                                                                                                           

    // Surface contact constraints
    const frictionParams FRICTION_PARAMS(MU, GAMMA, FOOT_WIDTH, FOOT_LENGTH);
    constraintsPtrVec_.push_back(std::make_unique<SurfaceContactConstraints>(cidPtr_, 
                                                                             FRICTION_PARAMS));
    constraintsNameVec_.push_back("contact constraints");

    // kinematics constraints
    constraintsPtrVec_.push_back(std::make_unique<DigitCustomizedConstraints>(model_input, 
                                                                              jtype_input, 
                                                                              trajPtr_, 
                                                                              dcidPtr_->ddcPtr_,
                                                                              gp_input));    
    constraintsNameVec_.push_back("customized constraints");            

    // periodic reset map constraints
    constraintsPtrVec_.push_back(std::make_unique<DigitSingleStepPeriodicityConstraints>(trajPtr_, 
                                                                                         cidPtr_,
                                                                                         FRICTION_PARAMS));    
    constraintsNameVec_.push_back("reset map constraints");     

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
bool DigitSingleStepOptimizer::get_nlp_info(
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

// [TNLP_eval_f]
// returns the value of the objective function
bool DigitSingleStepOptimizer::eval_f(
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

    cidPtr_->compute(z, false);

    obj_value = 0;

    // minimize control torque
    for ( Index i = 0; i < cidPtr_->N; i++ ) {
        obj_value += sqrt(cidPtr_->tau(i).dot(cidPtr_->tau(i)));
    }
    obj_value /= cidPtr_->N;

    // minimize initial velocity
    const VecX& initial_velocity = cidPtr_->trajPtr_->q_d(0);
    obj_value += 40 * sqrt(initial_velocity.dot(initial_velocity));

    // minimize initial acceleration
    const VecX& initial_acceleration = cidPtr_->trajPtr_->q_dd(0);
    obj_value += 20 * sqrt(initial_acceleration.dot(initial_acceleration));

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DigitSingleStepOptimizer::eval_grad_f(
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

    cidPtr_->compute(z, true);

    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] = 0;
    }

    for ( Index i = 0; i < cidPtr_->N; i++ ) {
        VecX v = cidPtr_->ptau_pz(i).transpose() * cidPtr_->tau(i);
        double norm = sqrt(cidPtr_->tau(i).dot(cidPtr_->tau(i)));   

        for ( Index j = 0; j < n; j++ ) {
            grad_f[j] += v(j) / norm;
        }
    }
    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] /= cidPtr_->N;
    }

    const VecX& initial_velocity = cidPtr_->trajPtr_->q_d(0);
    const VecX& initial_velocity_pz = cidPtr_->trajPtr_->pq_d_pz(0).transpose() * initial_velocity;
    const double initial_velocity_norm = sqrt(initial_velocity.dot(initial_velocity));
    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] += 40 * initial_velocity_pz(i) / initial_velocity_norm;
    }

    const VecX& initial_acceleration = cidPtr_->trajPtr_->q_dd(0);
    const VecX& initial_acceleration_pz = cidPtr_->trajPtr_->pq_dd_pz(0).transpose() * initial_acceleration;
    const double initial_acceleration_norm = sqrt(initial_acceleration.dot(initial_acceleration));
    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] += 20 * initial_acceleration_pz(i) / initial_acceleration_norm;
    }

    return true;
}
// [TNLP_eval_grad_f]

}; // namespace Digit
}; // namespace IDTO