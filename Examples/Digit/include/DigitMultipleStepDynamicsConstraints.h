
#ifndef DIGIT_MULTIPLE_STEP_CONSTRAINTS_H
#define DIGIT_MULTIPLE_STEP_CONSTRAINTS_H

#include "DigitDynamicsConstraints.h"

namespace IDTO {
namespace Digit {

// We always assume that the number of walking steps is even here!!!
class DigitMultipleStepDynamicsConstraints : public DigitDynamicsConstraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitMultipleStepDynamicsConstraints() = default;

    // Constructor
    DigitMultipleStepDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input,
                                         const Eigen::VectorXi& jtype_input, 
                                         char stanceLeg, 
                                         const Transform& stance_foot_T_des_input,
                                         const int N_input,
                                         const int NSteps_input);

    // Destructor
    ~DigitMultipleStepDynamicsConstraints() = default;

    // class methods:
    void setupJointPosition(VecX& q, bool compute_derivatives = true) final override;

    void setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives = true) final override;

    // class members:
        // number of time instances to evaluate constraints for each walking step
    int N = 0; 

        // number of walking steps in the gait optimization problem
    int NSteps = 0;

        // counter to record the number of times setupJointPosition is called
    int counter1 = 0;
    int counter2 = 0;
    
        // flag to indicate the first call of setupJointPosition
    bool firstCall = true;
};

int fillDependent_f(const gsl_vector* x, void *params, gsl_vector* f);

int fillDependent_df(const gsl_vector* x, void *params, gsl_matrix* J);

int fillDependent_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J);

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_MULTIPLE_STEP_CONSTRAINTS_H
