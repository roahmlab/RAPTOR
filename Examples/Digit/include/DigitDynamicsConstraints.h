
#ifndef DIGIT_CONSTRAINTS_H
#define DIGIT_CONSTRAINTS_H

#include <string>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "Utils.h"
#include "DigitConstants.h"
#include "DynamicsConstraints.h"
#include "ForwardKinematics.h"

namespace IDTO {
namespace Digit {


const std::string dependentJointNames[NUM_DEPENDENT_JOINTS] =
                                       {"Px",
                                        "Py",
                                        "Pz",
                                        "Rx",
                                        "Ry",
                                        "Rz",
                                        "left_achilles_rod",
                                        "left_ach2",
                                        "left_tarsus",
                                        "left_toe_A_rod",
                                        "left_A2",
                                        "left_toe_B_rod",
                                        "left_B2",
                                        "left_toe_pitch",
                                        "left_toe_roll",
                                        "right_achilles_rod",
                                        "right_ach2",
                                        "right_tarsus",
                                        "right_toe_A_rod",
                                        "right_A2",
                                        "right_toe_B_rod",
                                        "right_B2",
                                        "right_toe_pitch",
                                        "right_toe_roll"};
                                        
const std::string independentJointNames[NUM_INDEPENDENT_JOINTS] =
                                       {"left_hip_roll",
                                        "left_hip_yaw",
                                        "left_hip_pitch",
                                        "left_knee",
                                        "left_toe_A",
                                        "left_toe_B",
                                        "right_hip_roll",
                                        "right_hip_yaw",
                                        "right_hip_pitch",
                                        "right_knee",
                                        "right_toe_A",
                                        "right_toe_B"};

class DigitDynamicsConstraints : public DynamicsConstraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitDynamicsConstraints() = default;

    // Constructor
    DigitDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input,
                             const Eigen::VectorXi& jtype_input, 
                             char stanceLeg, 
                             const Transform& stance_foot_T_des_input);

    // Destructor
    ~DigitDynamicsConstraints() = default;

    // class methods:
        // swap the stance leg for reset map constraint evaluation
    void reinitialize() override;

        // return the index of id th dependent joint
    int return_dependent_joint_index(const int id) override;

        // return the index of id th independent joint
    int return_independent_joint_index(const int id) override;

        // fill in dependent indeces in a vector
    void fill_dependent_vector(VecX& r, const VecX& v, const bool setZero = false) override;

        // fill in independent indeces in a vector
    void fill_independent_vector(VecX& r, const VecX& v, const bool setZero = false) override;

        // fill in dependent columns in a matrix
    void fill_dependent_columns(MatX& r, const MatX& m, const bool setZero = false) override;

        // fill in independent columns in a matrix
    void fill_independent_columns(MatX& r, const MatX& m, const bool setZero = false) override;

        // fill in dependent rows in a matrix
    void fill_dependent_rows(MatX& r, const MatX& m, const bool setZero = false) override;

        // fill in independent rows in a matrix
    void fill_independent_rows(MatX& r, const MatX& m, const bool setZero = false) override;

        // return dependent indeces in a vector
    VecX get_dependent_vector(const VecX& v) override;

        // return independent indeces in a vector
    VecX get_independent_vector(const VecX& v) override;

        // return dependent columns in a matrix
    void get_dependent_columns(MatX& r, const MatX& m) override;

        // return independent columns in a matrix
    void get_independent_columns(MatX& r, const MatX& m) override;

        // return dependent rows in a matrix
    void get_dependent_rows(MatX& r, const MatX& m) override;

        // return independent rows in a matrix
    void get_independent_rows(MatX& r, const MatX& m) override;

        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    void setupJointPosition(VecX& q, bool compute_derivatives = true) override;

        // constraint c(q)
    void get_c(const VecX& q) override;

        // J = jacobian(c, q)
    void get_J(const VecX& q) override;
    
        // Jx_partial_dq = jacobian(J * x, q)
        // where x is a vector that is not dependent on q
    void get_Jx_partial_dq(const VecX& q, const VecX& x) override;

        // JTx_partial_dq = jacobian(J^T * x, q)
        // where x is a vector that is not dependent on q
    void get_JTx_partial_dq(const VecX& q, const VecX& x) override;

        // Jxy_partial_dq = jacobian(jacobian(J * x, q) * y, q)
        // where x and y are vectors that are not dependent on q
    void get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) override;

    // class members:
    std::shared_ptr<Model> modelPtr_ = nullptr;

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_ = nullptr;

        // dep/indep joint indeces
    int dependentJointIds[NUM_DEPENDENT_JOINTS] = {0};
    int independentJointIds[NUM_INDEPENDENT_JOINTS] = {0};

        // jtype copy
    Eigen::VectorXi jtype;

        // for contact constraints (forward kinematics mainly)
    Model::JointIndex contact_joint_id = 0;

        // forward kinematics derivatives
    Eigen::Array<MatX, 3, 1> H_translation;
    Eigen::Array<MatX, 3, 1> H_rotation;

    Eigen::Array<MatX, 3, 1> TOx_translation;
    Eigen::Array<MatX, 3, 1> TOx_rotation;

        // contact foot transforms
    char stanceLeg = 'L';

    Transform stance_foot_T;
    Transform stance_foot_endT;

    Transform stance_foot_T_des;

    VecX qcopy;
};

int fillDependent_f(const gsl_vector* x, void *params, gsl_vector* f);

int fillDependent_df(const gsl_vector* x, void *params, gsl_matrix* J);

int fillDependent_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J);

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_CONSTRAINTS_H
