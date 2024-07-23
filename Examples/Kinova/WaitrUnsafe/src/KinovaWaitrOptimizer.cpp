#include "KinovaWaitrOptimizer.h"

namespace RAPTOR {
namespace Kinova {

// // constructor
// KinovaWaitrOptimizer::KinovaWaitrOptimizer()
// {
// }


// // destructor
// KinovaWaitrOptimizer::~KinovaWaitrOptimizer()
// {
// }

// [TNLP_set_parameters]
bool KinovaWaitrOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    const WaitrTrajectoryParameters& atp_input,
    const contactSurfaceParams& csp_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input,
    const VecX& qdes_input,
    const int tplan_n_input,
    const VecX& joint_limits_buffer_input,
    const VecX& velocity_limits_buffer_input,
    const VecX& torque_limits_buffer_input
 ) 
{
    x0 = x0_input;
    qdes = qdes_input;
    tplan_n = tplan_n_input;

    trajPtr_ = std::make_shared<WaitrBezierCurves>(T_input, 
                                                   N_input, 
                                                   model_input.nq - 1, 
                                                   Chebyshev, 
                                                   atp_input);
                                                   
    idPtr_ = std::make_shared<CustomizedInverseDynamics>(model_input,
                                                         jtype_input,
                                                         trajPtr_);
    
    // read joint limits from KinovaConstants.h
    VecX JOINT_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS) + 
                                  joint_limits_buffer_input;

    VecX JOINT_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS) -
                                  joint_limits_buffer_input;

    // read velocity limits from KinovaConstants.h
    VecX VELOCITY_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(VELOCITY_LIMITS_LOWER, NUM_JOINTS) + 
                                     velocity_limits_buffer_input;

    VecX VELOCITY_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(VELOCITY_LIMITS_UPPER, NUM_JOINTS) -
                                     velocity_limits_buffer_input;

    // read torque limits from KinovaConstants.h
    VecX TORQUE_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_LOWER, NUM_JOINTS) + 
                                   torque_limits_buffer_input;

    VecX TORQUE_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(TORQUE_LIMITS_UPPER, NUM_JOINTS) -
                                   torque_limits_buffer_input;

    // Joint limits
    constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
                                                               JOINT_LIMITS_LOWER_VEC, 
                                                               JOINT_LIMITS_UPPER_VEC));
    constraintsNameVec_.push_back("joint limits");

    // Velocity limits
    constraintsPtrVec_.push_back(std::make_unique<VelocityLimits>(trajPtr_, 
                                                                  VELOCITY_LIMITS_LOWER_VEC, 
                                                                  VELOCITY_LIMITS_UPPER_VEC));
    constraintsNameVec_.push_back("velocity limits");        

    // Torque limits
    constraintsPtrVec_.push_back(std::make_unique<TorqueLimits>(trajPtr_, 
                                                                idPtr_,
                                                                TORQUE_LIMITS_LOWER_VEC, 
                                                                TORQUE_LIMITS_UPPER_VEC));
    constraintsNameVec_.push_back("torque limits");   

    // Contact constraints with the object
    constraintsPtrVec_.push_back(std::make_unique<WaitrContactConstraints>(idPtr_, 
                                                                           csp_input));                                                         
    constraintsNameVec_.push_back("contact constraints");

    // Customized constraints (collision avoidance with obstacles)
    constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_,
                                                                               model_input,
                                                                               jtype_input,
                                                                               boxCenters_input,
                                                                               boxOrientation_input,
                                                                               boxSize_input));   
    constraintsNameVec_.push_back("obstacle avoidance constraints");                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                                                                        
    assert(x0.size() == trajPtr_->varLength);
    assert(qdes.size() == trajPtr_->Nact);
    assert(tplan_n >= 0 && tplan_n < trajPtr_->N);

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
// returns some info about the nlp
bool KinovaWaitrOptimizer::get_nlp_info(
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
bool KinovaWaitrOptimizer::eval_f(
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

    const VecX& qplan = trajPtr_->q(tplan_n);

    obj_value = pow(Utils::wrapToPi(qplan[0] - qdes[0]), 2) + // These are continuous joints
                pow(Utils::wrapToPi(qplan[2] - qdes[2]), 2) + 
                pow(Utils::wrapToPi(qplan[4] - qdes[4]), 2) + 
                pow(Utils::wrapToPi(qplan[6] - qdes[6]), 2) + 
                pow(qplan[1] - qdes[1], 2) +           // These are not continuous joints
                pow(qplan[3] - qdes[3], 2) + 
                pow(qplan[5] - qdes[5], 2);

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool KinovaWaitrOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, true);

    const VecX& qplan = trajPtr_->q(tplan_n);
    const MatX& pqplan_pz = trajPtr_->pq_pz(tplan_n);

    for(Index i = 0; i < n; i++){
        if (i % 2 == 0) {
            grad_f[i] = (2 * Utils::wrapToPi(qplan[i] - qdes[i]) * pqplan_pz(i, i));
        }
        else {
            grad_f[i] = (2 * (qplan[i] - qdes[i]) * pqplan_pz(i, i));
        }
    }

    return true;
}
// [TNLP_eval_grad_f]

}; // namespace Kinova
}; // namespace RAPTOR