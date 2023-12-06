#include "DigitSingleStepOptimizer.h"

namespace IDTO {
namespace Digit {

using std::cout;
using std::endl;

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
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    const GaitParameters& gp_input
 ) 
{
    x0 = x0_input;

    fcPtr_ = std::make_shared<FourierCurves>(T_input, 
                                             N_input, 
                                             NUM_INDEPENDENT_JOINTS, 
                                             Chebyshev, 
                                             degree_input);

        // convert to their base class pointers
    trajPtr_ = fcPtr_;

    // stance foot is left foot by default
    char stanceLeg = 'L';
    Transform stance_foot_T_des(3, -M_PI / 2);
    dcidPtr_ = std::make_shared<DigitConstrainedInverseDynamics>(model_input, 
                                                                 trajPtr_,
                                                                 NUM_DEPENDENT_JOINTS, 
                                                                 jtype_input, 
                                                                 stanceLeg, 
                                                                 stance_foot_T_des);                                                          

    
    // convert joint limits from degree to radian
    VecX JOINT_LIMITS_LOWER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_LOWER_VEC(i) = deg2rad(JOINT_LIMITS_LOWER[i]);
    }

    // convert joint limits from degree to radian   
    VecX JOINT_LIMITS_UPPER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_UPPER_VEC(i) = deg2rad(JOINT_LIMITS_UPPER[i]);
    }

    VecX TORQUE_LIMITS_LOWER_VEC(NUM_INDEPENDENT_JOINTS);
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        TORQUE_LIMITS_LOWER_VEC(i) = TORQUE_LIMITS_LOWER[i];
    }

    VecX TORQUE_LIMITS_UPPER_VEC(NUM_INDEPENDENT_JOINTS);
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        TORQUE_LIMITS_UPPER_VEC(i) = TORQUE_LIMITS_UPPER[i];
    }

    // Joint limits
        // convert to their base class pointers
    dcPtr_ = dcidPtr_->dcPtr_;
    constraintsPtrVec_.push_back(std::make_unique<ConstrainedJointLimits>(trajPtr_, 
                                                                          dcPtr_, 
                                                                          JOINT_LIMITS_LOWER_VEC, 
                                                                          JOINT_LIMITS_UPPER_VEC));
    constraintsScale.push_back(1.0);                                                                      

    // Torque limits
        // convert to their base class pointers
    idPtr_ = dcidPtr_;
    constraintsPtrVec_.push_back(std::make_unique<TorqueLimits>(trajPtr_, 
                                                                idPtr_, 
                                                                TORQUE_LIMITS_LOWER_VEC, 
                                                                TORQUE_LIMITS_UPPER_VEC));  
    constraintsScale.push_back(1.0);                                                            

    // Surface contact constraints
        // convert to their base class pointers
    cidPtr_ = dcidPtr_;
    constraintsPtrVec_.push_back(std::make_unique<SurfaceContactConstraints>(cidPtr_, 
                                                                             MU, 
                                                                             GAMMA, 
                                                                             FOOT_WIDTH,
                                                                             FOOT_LENGTH));
    constraintsScale.push_back(1.0);    

    // Other constraints for gait optimization
        // convert to their base class pointers
    constraintsPtrVec_.push_back(std::make_unique<DigitCustomizedConstraints>(model_input, 
                                                                              jtype_input, 
                                                                              trajPtr_, 
                                                                              dcPtr_,
                                                                              gp_input));
    constraintsScale.push_back(1.0);                                                                                                                                                                                                                                                                                                                                                                                                       

    assert(x0.size() == trajPtr_->varLength);

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
       throw std::runtime_error("*** Error wrong value of n in eval_f!");
    }

    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    dcidPtr_->compute(z, false);

    obj_value = 0;
    for ( Index i = 0; i < dcidPtr_->N; i++ ) {
        obj_value += dcidPtr_->tau(i).dot(dcidPtr_->tau(i));
    }

    // obj_value *= 1e-3;

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
       throw std::runtime_error("*** Error wrong value of n in eval_f!");
    }

    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    dcidPtr_->compute(z, true);

    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] = 0;
    }

    for ( Index i = 0; i < dcidPtr_->N; i++ ) {
        VecX v = 2 * dcidPtr_->ptau_pz(i).transpose() * dcidPtr_->tau(i);

        for ( Index j = 0; j < n; j++ ) {
            grad_f[j] += v(j);
        }
    }

    return true;
}
// [TNLP_eval_grad_f]

}; // namespace Digit
}; // namespace IDTO