#ifndef ENDEFFECTOR_REGRESSOR_CONDITION_NUMBER
#define ENDEFFECTOR_REGRESSOR_CONDITION_NUMBER

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "RegressorInverseDynamics.h"
#include "Utils.h"

namespace RAPTOR {

namespace Torque {
constexpr double SQUARE_ROOT_THRESHOLD = 1e-6;
}; // namespace Torque

/*
minimize the condition number of the torque regressor matrix over the entire trajectory
*/
class EndEffectorRegressorConditionNumber : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    EndEffectorRegressorConditionNumber() = default;

    // Constructor
    EndEffectorRegressorConditionNumber(std::shared_ptr<Trajectories>& trajPtr_input, 
                                        std::shared_ptr<RegressorInverseDynamics>& ridPtr_input);

    // Destructor
    ~EndEffectorRegressorConditionNumber() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
};

}; // namespace RAPTOR

#endif // ENDEFFECTOR_REGRESSOR_CONDITION_NUMBER
