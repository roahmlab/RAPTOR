#include "DigitSingleStepOptimizer.h"

namespace RAPTOR {
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
    const GaitParameters& gp_input,
    const char stanceLeg,
    const Transform& stance_foot_T_des,
    bool periodic,
    const VecX q0_input,
    const VecX q_d0_input
) 
{
    x0 = x0_input;
                                          
    bcPtr_ = std::make_shared<BezierCurves>(T_input, 
                                            N_input, 
                                            NUM_INDEPENDENT_JOINTS, 
                                            time_discretization_input, 
                                            degree_input);     

    if (q0_input.size() == NUM_INDEPENDENT_JOINTS) {
        bcPtr_->constrainInitialPosition(q0_input);
    }
    if (q_d0_input.size() == NUM_INDEPENDENT_JOINTS) {
        bcPtr_->constrainInitialVelocity(q_d0_input);
    }

    // convert to base class
    trajPtr_ = bcPtr_;

    // add v_reset and lambda_reset to the end of the decision variables                                         
    trajPtr_->varLength += NUM_JOINTS + NUM_DEPENDENT_JOINTS;

    dcidPtr_ = std::make_shared<DigitConstrainedInverseDynamics>(model_input, 
                                                                 trajPtr_,
                                                                 NUM_DEPENDENT_JOINTS, 
                                                                 stanceLeg, 
                                                                 stance_foot_T_des);                                                          
    cidPtr_ = dcidPtr_; // convert to base class
    idPtr_ = cidPtr_; // convert to base class
    
    // convert joint limits from degree to radian
    VecX JOINT_LIMITS_LOWER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS));    
    VecX JOINT_LIMITS_UPPER_VEC = 
        Utils::deg2rad(
            Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS));

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
    const rectangleContactSurfaceParams FRICTION_PARAMS(MU, GAMMA, FOOT_WIDTH, FOOT_LENGTH);
    constraintsPtrVec_.push_back(std::make_unique<RectangleSurfaceContactConstraints>(cidPtr_, 
                                                                                      FRICTION_PARAMS));
    constraintsNameVec_.push_back("contact constraints");

    // kinematics constraints
    constraintsPtrVec_.push_back(std::make_unique<DigitCustomizedConstraints>(model_input, 
                                                                              trajPtr_, 
                                                                              dcidPtr_->ddcPtr_,
                                                                              gp_input));    
    constraintsNameVec_.push_back("customized constraints");            

    // periodic reset map constraints
    if (periodic) {
        constraintsPtrVec_.push_back(std::make_unique<DigitSingleStepPeriodicityConstraints>(trajPtr_, 
                                                                                             dcidPtr_,
                                                                                             FRICTION_PARAMS));    
        constraintsNameVec_.push_back("reset map constraints");     
    }

    // Cost functions
    costsPtrVec_.push_back(std::make_unique<MinimizePower>(trajPtr_, 
                                                           idPtr_));
    costsWeightVec_.push_back(1.0);
    costsNameVec_.push_back("minimize power");

    costsPtrVec_.push_back(std::make_unique<MinimizeInitialVelocity>(trajPtr_));
    costsWeightVec_.push_back(100.0);
    costsNameVec_.push_back("minimize initial velocity");

    costsPtrVec_.push_back(std::make_unique<MinimizeInitialAcceleration>(trajPtr_));
    costsWeightVec_.push_back(20.0);
    costsNameVec_.push_back("minimize initial acceleration");

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

}; // namespace Digit
}; // namespace RAPTOR