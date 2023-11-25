
#ifndef DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
#define DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

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

    // class methods:
        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    void setupJointPosition(VecX& q) override;
};

int fillDependent_f(const gsl_vector* x, void *params, gsl_vector* f);

int fillDependent_df(const gsl_vector* x, void *params, gsl_matrix* J);

int fillDependent_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J);

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
