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
    using MatX = Eigen::MatrixXd;
    using QRSolver = Eigen::ColPivHouseholderQR<MatX>;

    // Constructor
    ConstrainedInverseDynamics() = default;

    // Constructor
    ConstrainedInverseDynamics(const Model& model_input, 
                               int N_input, 
                               int numDependentJoints_input,
                               std::unique_ptr<DynamicsConstraints>& dynamics_constraints_input);

    // Destructor
    ~ConstrainedInverseDynamics() = default;

    // class methods:
        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    virtual void setupJointPosition(VecX& q) = 0;

        // fill in dependent joint positions and velocities in the full joint vector q and v
        // that satisfies the constraints
    // virtual void setupJointPositionVelocity(VecX& q, VecX& v);

        // fill in dependent joint positions, velocities, and accelerations in the full joint vector q, v, and a
        // that satisfies the constraints
    virtual void setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives = true);
    
    virtual void compute(Eigen::Array<VecX, 1, Eigen::Dynamic>& q, 
                         Eigen::Array<VecX, 1, Eigen::Dynamic>& v, 
                         Eigen::Array<VecX, 1, Eigen::Dynamic>& a,
                         bool compute_derivatives = true);

    // class members:
    int numDependentJoints = 0;
    int numIndependentJoints = 0;

        // declare a DynamicsConstraints instance outside of this class
        // use a shared pointer here to avoid copying
    std::shared_ptr<DynamicsConstraints> dynamicsConstraintsPtr_;

        // updated in setupJointPositionVelocityAcceleration()
    QRSolver J_dep_qr;
    QRSolver J_dep_T_qr;

        // updated in setupJointPositionVelocityAcceleration()
    MatX J_dep;
    MatX J_indep;

        // updated in compute()
    VecX tau_dep;
    VecX tau_indep;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> lambda;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pq;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pv;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pa;
};

}; // namespace IDTO

#endif // CONSTRAINEDINVERSEDYNAMICS_H
