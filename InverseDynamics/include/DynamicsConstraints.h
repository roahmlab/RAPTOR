#ifndef DYNAMICSCONSTRAINTS_H
#define DYNAMICSCONSTRAINTS_H

#include "InverseDynamics.h"

namespace IDTO {

class DynamicsConstraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using QRSolver = Eigen::ColPivHouseholderQR<MatX>;

    // Constructor
    DynamicsConstraints() = default;

    // Constructor
    DynamicsConstraints(const Model& model_input, int numDependentJoints_input);

    // Destructor
    ~DynamicsConstraints() = default;

    // class methods:
    virtual void reinitialize();

    // You need to implement all following methods in your derived class!!!
        // return the index of id th dependent joint
    virtual int return_dependent_joint_index(const int id) = 0;

        // return the index of id th independent joint
    virtual int return_independent_joint_index(const int id) = 0;

        // fill in dependent indeces in a vector
    virtual void fill_dependent_vector(VecX& r, const VecX& v, const bool setZero = false) = 0;

        // fill in independent indeces in a vector
    virtual void fill_independent_vector(VecX& r, const VecX& v, const bool setZero = false) = 0;

        // fill in dependent columns in a matrix
    virtual void fill_dependent_columns(MatX& r, const MatX& m, const bool setZero = false) = 0;

        // fill in independent columns in a matrix
    virtual void fill_independent_columns(MatX& r, const MatX& m, const bool setZero = false) = 0;

        // fill in dependent rows in a matrix
    virtual void fill_dependent_rows(MatX& r, const MatX& m, const bool setZero = false) = 0;

        // fill in independent rows in a matrix
    virtual void fill_independent_rows(MatX& r, const MatX& m, const bool setZero = false) = 0;

        // return dependent indeces in a vector
    virtual VecX get_dependent_vector(const VecX& v) = 0;

        // return independent indeces in a vector
    virtual VecX get_independent_vector(const VecX& v) = 0;

        // return dependent columns in a matrix
    virtual void get_dependent_columns(MatX& r, const MatX& m) = 0;

        // return independent columns in a matrix
    virtual void get_independent_columns(MatX& r, const MatX& m) = 0;

        // return dependent rows in a matrix
    virtual void get_dependent_rows(MatX& r, const MatX& m) = 0;

        // return independent rows in a matrix
    virtual void get_independent_rows(MatX& r, const MatX& m) = 0;

    // class methods:
        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    virtual void setupJointPosition(VecX& q, bool compute_derivatives = true) = 0;

        // fill in dependent joint positions and velocities in the full joint vector q and v
        // that satisfies the constraints
    // virtual void setupJointPositionVelocity(VecX& q, VecX& v);

        // fill in dependent joint positions, velocities, and accelerations in the full joint vector q, v, and a
        // that satisfies the constraints
    virtual void setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives = true);

    // constraint c(q)
    virtual void get_c(const VecX& q) = 0;

    // J = jacobian(c, q)
    virtual void get_J(const VecX& q) = 0;
    
    // Jx_partial_dq = jacobian(J * x, q)
    // where x is a vector that is not dependent on q
    virtual void get_Jx_partial_dq(const VecX& q, const VecX& x) = 0;

    // JTx_partial_dq = jacobian(J^T * x, q)
    // where x is a vector that is not dependent on q
    virtual void get_JTx_partial_dq(const VecX& q, const VecX& x) = 0;

    // Jxy_partial_dq = jacobian(jacobian(J * x, q) * y, q)
    // where x and y are vectors that are not dependent on q
    virtual void get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) = 0;

    // class members:
    std::unique_ptr<Model> modelPtr_;

    int numDependentJoints = 0;
    int numIndependentJoints = 0;

        // compute results are stored here
    VecX c;
    MatX J;
    MatX Jx_partial_dq;
    MatX JTx_partial_dq;
    MatX Jxy_partial_dq;

    MatX pq_dep_pq_indep;
    MatX pv_dep_pq;
    MatX pv_dep_pv_indep;
    MatX pa_dep_pq;
    MatX pa_dep_pv;
    MatX pa_dep_pa_indep;

        // updated in setupJointPositionVelocityAcceleration()
    QRSolver J_dep_qr;
    QRSolver J_dep_T_qr;

        // updated in setupJointPositionVelocityAcceleration()
    MatX J_dep;
    MatX J_indep;

        // temporary variables updated in setupJointPositionVelocityAcceleration()
    MatX P_dep;
    VecX Pa_indep;
    MatX temp1;
    VecX temp2_1;
    MatX temp2;
    MatX temp3;
};

}; // namespace IDTO

#endif // DYNAMICSCONSTRAINTS_H
