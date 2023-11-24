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

    // Constructor
    DynamicsConstraints() = default;

    // Constructor
    DynamicsConstraints(const Model& model_input, int numIndependentJoints);

    // Destructor
    ~DynamicsConstraints() = default;

    // class methods:
    // You need to implement all following methods in your derived class!!!

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

        // compute results are stored here
    VecX c;
    MatX J;
    MatX Jx_partial_dq;
    MatX JTx_partial_dq;
    MatX Jxy_partial_dq;
};

}; // namespace IDTO

#endif // DYNAMICSCONSTRAINTS_H
