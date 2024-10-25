#ifndef DIGIT_PARAMETERS_IDENTIFICATION_H
#define DIGIT_PARAMETERS_IDENTIFICATION_H

#include "Optimizer.h"
#include "LMIConstraints.h"
#include "LinearConstraints.h"
#include "DigitWholeBodyDynamicsConstraints.h"

#include "pinocchio/algorithm/regressor.hpp"
#include "pinocchio/algorithm/center-of-mass.hpp"

namespace RAPTOR {
namespace DigitWholeBodySysID {

class DigitSystemIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Mat3 = Eigen::Matrix3d;

    /** Default constructor */
    DigitSystemIdentification() = default;

    /** Default destructor */
    ~DigitSystemIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const std::shared_ptr<MatX>& posPtr_input,
        const std::shared_ptr<MatX>& torquePtr_input);

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
    bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    ) final override;

    /** Method to return the bounds for my problem */
    bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
    ) final override;

    /** Method to return the objective value */
    bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    ) final override;

    /** Method to return the gradient of the objective */
    bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    ) final override;

    /** Method to return the hessian of the objective */
    bool eval_hess_f(
        Index         n,
        const Number* x,
        bool          new_x,
        MatX&         hess_f
    ) final override;
    
    /**@name Methods to block default compiler methods.
    *
    * The compiler automatically generates the following three methods.
    *  Since the default compiler implementation is generally not what
    *  you want (for all but the most simple classes), we usually
    *  put the declarations of these methods in the private section
    *  and never implement them. This prevents the compiler from
    *  implementing an incorrect "default" behavior without us
    *  knowing. (See Scott Meyers book, "Effective C++")
    */
    //@{
    DigitSystemIdentification(
       const DigitSystemIdentification&
    );

    DigitSystemIdentification& operator=(
       const DigitSystemIdentification&
    );

    const double default_maximum_uncertainty = 0.5; // default maximum uncertainty

    std::shared_ptr<Model> modelPtr_; // robot model
    std::shared_ptr<Data> dataPtr_; // robot data

    std::shared_ptr<DynamicsConstraints> ddcPtr_; // digit dynamics constraints

    std::vector<int> nontrivialLinkIds; // nontrivial link ids

    MatX A;
    VecX b;

    VecX tau_fixed; // torque corresponding to link inertial parameters that are fixed for the optimization

    MatX FullObservationMatrix; // full observation matrix
    VecX tau_estimated; // estimated torque

    Eigen::Array<MatX, 1, Eigen::Dynamic> Js; // jacobian matrices

    // shared pointers to data
    std::shared_ptr<MatX> posPtr_;
    std::shared_ptr<MatX> torquePtr_;

    MatX tau_inertials; // computed from the trajectory without friction

    int Nact = 0; // number of actuated joints
    int N = 0; // number of samples
};

class CustomizedConstraints : public Constraints {
public: 
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    CustomizedConstraints() = default;

    // Constructor
    CustomizedConstraints(const double totalMass_input,
                          const int N_input,
                          const int lambdaOffset_input,
                          const int varLength) :
        totalMass(totalMass_input),
        N(N_input),
        lambdaOffset(lambdaOffset_input),
        Constraints(N_input * 3, varLength) {};

    // Destructor
    ~CustomizedConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override;

        // print violation information
    virtual void print_violation_info() override;

    // class members:
    double totalMass = 0;
    int N = 0;
    int lambdaOffset = 0;
};

}; // namespace DigitWholeBodySysID
}; // namespace RAPTOR

#endif // DIGIT_PARAMETERS_IDENTIFICATION_H