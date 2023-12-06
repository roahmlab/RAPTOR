#include "KinovaOptimizer.h"

namespace IDTO {
namespace Kinova {

// // constructor
// KinovaOptimizer::KinovaOptimizer()
// {
// }


// // destructor
// KinovaOptimizer::~KinovaOptimizer()
// {
// }


bool KinovaOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input
 ) 
{
    x0 = x0_input;

    fcPtr_ = std::make_shared<FourierCurves>(T_input, 
                                             N_input, 
                                             model_input.nq, 
                                             Chebyshev, 
                                             degree_input);

        // convert to their base class pointers
    trajPtr_ = fcPtr_;
    
    // // convert joint limits from degree to radian
    // VecX JOINT_LIMITS_LOWER_VEC(NUM_JOINTS);
    // for (int i = 0; i < NUM_JOINTS; i++) {
    //     JOINT_LIMITS_LOWER_VEC(i) = deg2rad(JOINT_LIMITS_LOWER[i]);
    // }

    // // convert joint limits from degree to radian   
    // VecX JOINT_LIMITS_UPPER_VEC(NUM_JOINTS);
    // for (int i = 0; i < NUM_JOINTS; i++) {
    //     JOINT_LIMITS_UPPER_VEC(i) = deg2rad(JOINT_LIMITS_UPPER[i]);
    // }

    // // Joint limits
    //     // convert to their base class pointers
    // constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
    //                                                            JOINT_LIMITS_LOWER_VEC, 
    //                                                            JOINT_LIMITS_UPPER_VEC));
    // constraintsScale.push_back(1.0);     

    constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_));      

    Transform endT;
    VecX desiredEndEffectorPos(6);
    desiredEndEffectorPos << -0.6507, 0.1673, 0.7073, -2.5521, 1.1871, -1.9560;

    constraintsPtrVec_.push_back(std::make_unique<EndEffectorConstraints>(model_input,
                                                                          jtype_input,
                                                                          endT,
                                                                          "joint_7",
                                                                          trajPtr_,
                                                                          desiredEndEffectorPos));                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                        
    assert(x0.size() == trajPtr_->varLength);

    return true;
}
// [TNLP_set_parameters]

bool KinovaOptimizer::get_nlp_info(
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
bool KinovaOptimizer::eval_f(
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

    trajPtr_->compute(z, false);

    obj_value = trajPtr_->q_dd(trajPtr_->N / 2).dot(trajPtr_->q_dd(trajPtr_->N / 2));

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool KinovaOptimizer::eval_grad_f(
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

    trajPtr_->compute(z, true);

    VecX temp = 2 * trajPtr_->pq_dd_pz(trajPtr_->N / 2).transpose() * trajPtr_->q_dd(trajPtr_->N / 2);

    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] = temp(i);
    }

    return true;
}
// [TNLP_eval_grad_f]

}; // namespace Kinova
}; // namespace IDTO