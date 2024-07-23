
#ifndef DIGIT_MODIFIED_CONSTRAINTS_H
#define DIGIT_MODIFIED_CONSTRAINTS_H

#include <string>
#include <vector>

#include "DigitModifiedConstants.h"
#include "DynamicsConstraints.h"
#include "ForwardKinematics.h"

namespace RAPTOR {
namespace DigitModified {

const std::string dependentJointNames[NUM_DEPENDENT_JOINTS] =
                                       {"Px",
                                        "Py",
                                        "Pz",
                                        "Rx",
                                        "Ry",
                                        "Rz"};
                                        
const std::string independentJointNames[NUM_INDEPENDENT_JOINTS] =
                                       {"left_hip_roll",
                                        "left_hip_yaw",
                                        "left_hip_pitch",
                                        "left_knee",
                                        "left_tarsus",
                                        "left_toe_pitch",
                                        "left_toe_roll",
                                        "right_hip_roll",
                                        "right_hip_yaw",
                                        "right_hip_pitch",
                                        "right_knee",
                                        "right_tarsus",
                                        "right_toe_pitch",
                                        "right_toe_roll"};

class DigitModifiedDynamicsConstraints : public DynamicsConstraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitModifiedDynamicsConstraints() = default;

    // Constructor
    DigitModifiedDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input,
                                     const Eigen::VectorXi& jtype_input, 
                                     char stanceLeg, 
                                     const Transform& stance_foot_T_des_input);

    // Destructor
    ~DigitModifiedDynamicsConstraints() = default;

    // class methods:
        // swap the stance leg for reset map constraint evaluation
    void reinitialize();

        // return the index of id th dependent joint
    int return_dependent_joint_index(const int id) final override;

        // return the index of id th independent joint
    int return_independent_joint_index(const int id) final override;

        // fill in dependent indeces in a vector
    void fill_dependent_vector(VecX& r, const VecX& v, const bool setZero = false) final override;

        // fill in independent indeces in a vector
    void fill_independent_vector(VecX& r, const VecX& v, const bool setZero = false) final override;

        // fill in dependent columns in a matrix
    void fill_dependent_columns(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in independent columns in a matrix
    void fill_independent_columns(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in dependent rows in a matrix
    void fill_dependent_rows(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in independent rows in a matrix
    void fill_independent_rows(MatX& r, const MatX& m, const bool setZero = false) final override;

        // return dependent indeces in a vector
    VecX get_dependent_vector(const VecX& v) final override;

        // return independent indeces in a vector
    VecX get_independent_vector(const VecX& v) final override;

        // return dependent columns in a matrix
    void get_dependent_columns(MatX& r, const MatX& m) final override;

        // return independent columns in a matrix
    void get_independent_columns(MatX& r, const MatX& m) final override;

        // return dependent rows in a matrix
    void get_dependent_rows(MatX& r, const MatX& m) final override;

        // return independent rows in a matrix
    void get_independent_rows(MatX& r, const MatX& m) final override;

        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    void setupJointPosition(VecX& q, bool compute_derivatives = true) final override;

        // constraint c(q)
    void get_c(const VecX& q) final override;

        // J = jacobian(c, q)
    void get_J(const VecX& q) final override;
    
        // Jx_partial_dq = jacobian(J * x, q)
        // where x is a vector that is not dependent on q
    void get_Jx_partial_dq(const VecX& q, const VecX& x) final override;

        // JTx_partial_dq = jacobian(J^T * x, q)
        // where x is a vector that is not dependent on q
    void get_JTx_partial_dq(const VecX& q, const VecX& x) final override;

        // Jxy_partial_dq = jacobian(jacobian(J * x, q) * y, q)
        // where x and y are vectors that are not dependent on q
    void get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) final override;

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

        // contact foot transforms
    char stanceLeg = 'L';
    
    Transform stance_foot_endT;
    Transform stance_foot_T_des;

        // intermediate results
    Vec3 c_translation;
    Vec3 c_rotation;

    MatX J_translation;
    MatX J_rotation;

    Eigen::Array<MatX, 3, 1> H_translation;
    Eigen::Array<MatX, 3, 1> H_rotation;

    Eigen::Array<MatX, 3, 1> TOx_translation;
    Eigen::Array<MatX, 3, 1> TOx_rotation;
};

}; // namespace DigitModified
}; // namespace RAPTOR

#endif // DIGIT_MODIFIED_CONSTRAINTS_H
