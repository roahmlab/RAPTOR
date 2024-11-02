#include "EndeffectorConditionNumberOptimizer.h"



namespace RAPTOR {
namespace Kinova {

// // constructor
// EndeffectorConditionNumberOptimizer::EndeffectorConditionNumberOptimizer()
// {
// }


// // destructor
// EndeffectorConditionNumberOptimizer::~EndeffectorConditionNumberOptimizer()
// {
// }

bool EndeffectorConditionNumberOptimizer::set_parameters(
    const VecX& x0_input,
    const Number T_input,
    const int N_input,
    const int degree_input,
    const double base_frequency_input,
    const Model& model_input, 
    const VecX& joint_limits_buffer_input,
    const VecX& velocity_limits_buffer_input,
    const VecX& torque_limits_buffer_input,
    const bool include_gripper_or_not,
    const double collison_buffer_input,
    Eigen::VectorXi jtype_input
) {
    enable_hessian = false;

    // fixed frequency fourier curves with 0 initial velocity
    // trajPtr_ = std::make_shared<FixedFrequencyFourierCurves>(T_input, 
    //                                                          N_input, 
    //                                                          model_input.nv, 
    //                                                          TimeDiscretization::Uniform, 
    //                                                          degree_input,
    //                                                          base_frequency_input);

    // fixed frequency fourier curves with 0 initial velocity
    VecX q0_input(model_input.nv);
    q0_input << 1.001089876408351, 0.09140272042061115,  -1.648806446891836,   2.381092213417765,
                     1.822374826812066,  0.1466609489107418,  0.9315315991321746;
    VecX q_d0_input = VecX::Zero(model_input.nv);

    trajPtr_ = std::make_shared<FixedFrequencyFourierCurves>(T_input, 
                                                             N_input, 
                                                             model_input.nv, 
                                                             TimeDiscretization::Uniform, 
                                                             degree_input,
                                                             base_frequency_input,
                                                             q0_input,
                                                             q_d0_input);

    ridPtr_ = std::make_shared<RegressorInverseDynamics>(model_input, 
                                                         trajPtr_,
                                                         jtype_input);

    // add end effector regressor condition number into cost
    costsPtrVec_.push_back(std::make_unique<EndEffectorRegressorConditionNumber>(trajPtr_, 
                                                                                 ridPtr_));
    costsWeightVec_.push_back(1.0);                                                                             
    costsNameVec_.push_back("end effector regressor condition number");

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
                                                                ridPtr_,
                                                                TORQUE_LIMITS_LOWER_VEC, 
                                                                TORQUE_LIMITS_UPPER_VEC));
    constraintsNameVec_.push_back("torque limits"); 

    // Customized constraints (collision avoidance with ground)
    std::vector<Vec3> Center = {Vec3(0.0, 0.0, 0.15),
                                 Vec3(0.53, 0.49, 0.56),  // back wall
                                 Vec3(-0.39, -0.84, 0.56), // bar near the control
                                 Vec3(-0.39, -0.17, 0.56), //bar bewteen 10 and 20 change to wall
                                 Vec3(0.0, 0.0, 1.12), //ceiling
                                 Vec3(0.47, -0.09, 1.04) // top camera  
                                };    
    std::vector<Vec3> Orientation = {Vec3(0.0, 0.0, 0.0),
                                     Vec3(0.0, 0.0, 0.0),
                                     Vec3(0.0, 0.0, 0.0),
                                     Vec3(0.0, 0.0, 0.0),
                                     Vec3(0.0, 0.0, 0.0),
                                     Vec3(0.0, 0.0, 0.0)
                                    };
    std::vector<Vec3> Size = {Vec3(5.0, 5.0, 0.01),
                              Vec3(5.0, 0.05, 1.12),
                              Vec3( 0.05, 0.05, 1.12),
                              Vec3( 0.05, 1.28, 1.28),
                              Vec3( 5, 5, 0.05),
                              Vec3( 0.15, 0.15, 0.15)
                            };    
    constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_,
                                                                               model_input,
                                                                               Center,
                                                                               Orientation,
                                                                               Size,
                                                                               include_gripper_or_not,
                                                                               collison_buffer_input,
                                                                               jtype_input));  
    constraintsNameVec_.push_back("obstacle avoidance constraints"); 

    x0 = x0_input.head(trajPtr_->varLength);

    return true;
}

bool EndeffectorConditionNumberOptimizer::get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
) {
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

}; // namespace Kinova
}; // namespace RAPTOR