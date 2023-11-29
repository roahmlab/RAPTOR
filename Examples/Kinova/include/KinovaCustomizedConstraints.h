#ifndef KINOVA_CUSTOMIZED_CONSTRAINTS_H
#define KINOVA_CUSTOMIZED_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"

#include <memory>

namespace IDTO {
namespace Kinova {

class KinovaCustomizedConstraints : public Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    KinovaCustomizedConstraints() = default;

    // Constructor
    KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~KinovaCustomizedConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<Trajectories>& trajPtr_;
};

}; // namespace Kinova
}; // namespace IDTO

#endif // KINOVA_CUSTOMIZED_CONSTRAINTS_H
