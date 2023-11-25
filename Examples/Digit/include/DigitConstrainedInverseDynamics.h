
#ifndef DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
#define DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"

namespace IDTO {
namespace Digit {

class DigitConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    // Constructor
    DigitConstrainedInverseDynamics() = default;

    // Constructor
    DigitConstrainedInverseDynamics(const Model& model_input, 
                                    int N_input, 
                                    int numDependentJoints_input,
                                    std::unique_ptr<DynamicsConstraints>& dynamics_constraints_input);

    // Destructor
    ~DigitConstrainedInverseDynamics() = default;
};

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
