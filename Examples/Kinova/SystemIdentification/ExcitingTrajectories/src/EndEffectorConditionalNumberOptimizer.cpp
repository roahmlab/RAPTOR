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
        // const std::string& regroupMatrixFileName,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input,
        Eigen::VectorXi jtype_input
) {
    enable_hessian = false;
    x0 = x0_input;
    joint_num =model_input.nv;

    // fixed frequency fourier curves with 0 initial velocity
    trajPtr_ = std::make_shared<FixedFrequencyFourierCurves>(T_input, 
                                                             N_input, 
                                                             model_input.nv, 
                                                             TimeDiscretization::Uniform, 
                                                             degree_input,
                                                             base_frequency_input);

    ridPtr_ = std::make_shared<RegressorInverseDynamics>(model_input, 
                                                         trajPtr_,
                                                         jtype_input);

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

    // // Customized constraints (collision avoidance with ground)
    std::vector<Vec3> groundCenter = {Vec3(0.0, 0.0, 0.04)};
    std::vector<Vec3> groundOrientation = {Vec3(0.0, 0.0, 0.0)};
    std::vector<Vec3> groundSize = {Vec3(5.0, 5.0, 0.01)};
    constraintsPtrVec_.push_back(std::make_unique<KinovaCustomizedConstraints>(trajPtr_,
                                                                               model_input,
                                                                               groundCenter,
                                                                               groundOrientation,
                                                                               groundSize,
                                                                               0.0,
                                                                               jtype_input));   
    constraintsNameVec_.push_back("obstacle avoidance constraints"); 

    // // check dimensions of regroupMatrix
    // if (ridPtr_->Y.cols() != regroupMatrix.rows()) {
    //     throw std::invalid_argument("EndeffectorConditionNumberOptimizer: regroupMatrix has wrong dimensions!");
    // }

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

bool EndeffectorConditionNumberOptimizer::eval_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number&       obj_value
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    ridPtr_->compute(z, false);

    MatX ObservationMatrix = ridPtr_->Y;
    int numParams = ridPtr_->phi.size();

    // last numparams cols and row every numParams/joint_num rows  
    MatX lastCols = ObservationMatrix.rightCols(numParams/joint_num);
    MatX EndeffectorY =  lastCols(Eigen::seq(0, Eigen::last, joint_num-1), Eigen::all);
  
    Eigen::JacobiSVD<MatX> svd(EndeffectorY, 
                               Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VecX& singularValues = svd.singularValues();
    const MatX& U = svd.matrixU();
    const MatX& V = svd.matrixV();

    Number sigmaMax = singularValues(0);
    Number sigmaMin = singularValues(singularValues.size() - 1);

    // log of 2-norm condition number (sigmaMax / sigmaMin)
    obj_value = std::log(sigmaMax) - std::log(sigmaMin);

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool EndeffectorConditionNumberOptimizer::eval_grad_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number*       grad_f
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    ridPtr_->compute(z, true);

 
    MatX ObservationMatrix = ridPtr_->Y;
    int numParams = ridPtr_->phi.size();
    std::cout << numParams << "parems" <<std ::endl;
    
    // last numparams cols and row every numParams/joint_num rows  
    MatX lastCols = ObservationMatrix.rightCols(numParams/joint_num);
    MatX EndeffectorY =  lastCols(Eigen::seq(0, Eigen::last,  joint_num-1), Eigen::all);


    Eigen::JacobiSVD<MatX> svd(EndeffectorY, 
                               Eigen::ComputeThinU | Eigen::ComputeThinV);

    const VecX& singularValues = svd.singularValues();
    const MatX& U = svd.matrixU();
    const MatX& V = svd.matrixV();

    Index lastRow = singularValues.size() - 1;
    Number sigmaMax = singularValues(0);
    Number sigmaMin = singularValues(lastRow);

    // refer to https://j-towns.github.io/papers/svd-derivative.pdf
    // for analytical gradient of singular values
    for (Index i = 0; i < n; i++) {

        MatX gradObservationMatrix = ridPtr_->pY_pz(i);

        MatX lastCols = gradObservationMatrix.rightCols(numParams/joint_num);
        MatX gradEndeffectorY =  lastCols(Eigen::seq(0, Eigen::last,  joint_num -1), Eigen::all);


        Number gradSigmaMax = U.col(0).transpose()       * gradEndeffectorY * V.col(0);
        Number gradSigmaMin = U.col(lastRow).transpose() * gradEndeffectorY * V.col(lastRow);

        grad_f[i] = gradSigmaMax / sigmaMax - gradSigmaMin / sigmaMin;
    }

    return true;
}

}; // namespace Kinova
}; // namespace RAPTOR