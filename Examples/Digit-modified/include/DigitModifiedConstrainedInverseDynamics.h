
#ifndef DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H
#define DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "DigitModifiedDynamicsConstraints.h"

namespace IDTO {
namespace DigitModified {

class DigitModifiedConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    DigitModifiedConstrainedInverseDynamics() = default;

    // Constructor
    DigitModifiedConstrainedInverseDynamics(const Model& model_input, 
                                            std::shared_ptr<Trajectories>& trajPtr_input,
                                            int numDependentJoints_input,
                                            const Eigen::VectorXi& jtype_input, 
                                            char stanceLeg, 
                                            const Transform& stance_foot_T_des_input);

    // Destructor
    ~DigitModifiedConstrainedInverseDynamics() = default;
};

}; // namespace DigitModified
}; // namespace IDTO

#endif // DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H
