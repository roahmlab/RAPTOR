#ifndef CONSTRAINEDINVERSEDYNAMICS_H
#define CONSTRAINEDINVERSEDYNAMICS_H

#include "DynamicsConstraints.h"

namespace IDTO {

// Class declaration
class ConstrainedInverseDynamics : public InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = MatX;

    // Constructor
    ConstrainedInverseDynamics() = default;

    // Constructor
    ConstrainedInverseDynamics(const Model& model_input, 
                               std::shared_ptr<Trajectories>& trajPtr_input,
                               int N_input, 
                               int numDependentJoints_input);

    // Destructor
    ~ConstrainedInverseDynamics() = default;
    
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true);

    // class members:
    int numDependentJoints = 0;
    int numIndependentJoints = 0;

        // declare a DynamicsConstraints instance outside of this class
        // use a shared pointer here to avoid copying
    std::shared_ptr<DynamicsConstraints> dcPtr_;

        // updated in compute()
    VecX tau_dep;
    VecX tau_indep;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> lambda;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pz;

        // temporary variables updated in compute()
    MatX pq_pz;
    MatX pv_pz;
    MatX pa_pz;

    MatX prnea_pq_dep;
    MatX prnea_pq_indep;
    MatX prnea_pv_dep;
    MatX prnea_pv_indep;
    MatX prnea_pa_dep;
    MatX prnea_pa_indep;

    MatX plambda_pq;
    MatX plambda_pv;
    MatX plambda_pa;

    MatX JTlambda_partial_dq;
    MatX JTlambda_partial_dq_dep;
    MatX JTlambda_partial_dq_indep;
};

}; // namespace IDTO

#endif // CONSTRAINEDINVERSEDYNAMICS_H
