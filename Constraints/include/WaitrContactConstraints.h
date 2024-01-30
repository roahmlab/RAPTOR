#ifndef WAITR_CONTACT_CONSTRAINTS_H
#define WAITR_CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "InverseDynamicsWithFixedJoints.h"
#include "Utils.h"

namespace IDTO {

typedef struct contactSurfaceParams_  {
    // TODO: define the contact parameters here, like friction coefficient, etc.
    // double mu = 0.7;
    // double gamma = 0.7;
    // double Lx = 0;
    // double Ly = 0;

    // contactSurfaceParams_(double mu_input, 
    //                 double gamma_input, 
    //                 double Lx_input, 
    //                 double Ly_input) :
    //     mu(mu_input), 
    //     gamma(gamma_input), 
    //     Lx(Lx_input), 
    //     Ly(Ly_input) {}
} contactSurfaceParams;

class WaitrContactConstraints : public Constraints {
public:
    using Force = pinocchio::Data::Force;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    WaitrContactConstraints() = default;

    // Constructor
    WaitrContactConstraints(std::shared_ptr<InverseDynamicsWithFixedJoints>& idPtr_input,
                            const contactSurfaceParams& csp_input);

    // Destructor
    ~WaitrContactConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<InverseDynamicsWithFixedJoints> idPtr_;

    contactSurfaceParams fp;
};

}; // namespace IDTO

#endif // WAITR_CONTACT_CONSTRAINTS_H
