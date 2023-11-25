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
                               int N_input, 
                               int numDependentJoints_input);

    // Destructor
    ~ConstrainedInverseDynamics() = default;
    
    virtual void compute(Eigen::Array<VecX, 1, Eigen::Dynamic>& q, 
                         Eigen::Array<VecX, 1, Eigen::Dynamic>& v, 
                         Eigen::Array<VecX, 1, Eigen::Dynamic>& a,
                         bool compute_derivatives = true);

    // class members:
    int numDependentJoints = 0;
    int numIndependentJoints = 0;

        // declare a DynamicsConstraints instance outside of this class
        // use a shared pointer here to avoid copying
    std::unique_ptr<DynamicsConstraints> dynamicsConstraintsPtr_;

        // updated in compute()
    VecX tau_dep;
    VecX tau_indep;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> lambda;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pq;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pv;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pa;

        // temporary variables updated in compute()
    MatX rnea_partial_dq_dep;
    MatX rnea_partial_dq_indep;
    MatX rnea_partial_dv_dep;
    MatX rnea_partial_dv_indep;
    MatX rnea_partial_da_dep;
    MatX rnea_partial_da_indep;

    MatX JTlambda_partial_dq;
    MatX JTlambda_partial_dq_dep;
    MatX JTlambda_partial_dq_indep;
};

}; // namespace IDTO

#endif // CONSTRAINEDINVERSEDYNAMICS_H
