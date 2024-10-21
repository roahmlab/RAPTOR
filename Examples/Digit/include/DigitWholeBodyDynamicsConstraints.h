
#ifndef DIGIT_WHOLEBODY_CONSTRAINTS_H
#define DIGIT_WHOLEBODY_CONSTRAINTS_H

#include <string>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "Utils.h"
#include "DynamicsConstraints.h"

#include "pinocchio/algorithm/frames.hpp"

namespace RAPTOR {
namespace DigitWholeBodySysID {

constexpr int NUM_DEPENDENT_JOINTS = 24;
constexpr int NUM_INDEPENDENT_JOINTS = 20;

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
                                        "right_toe_B",
                                        "left_shoulder_roll",
                                        "left_shoulder_pitch",
                                        "left_shoulder_yaw",
                                        "left_elbow",
                                        "right_shoulder_roll",
                                        "right_shoulder_pitch",
                                        "right_shoulder_yaw",
                                        "right_elbow"};

class DigitWholeBodyDynamicsConstraints : public DynamicsConstraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    DigitWholeBodyDynamicsConstraints() = default;

    // Constructor
    DigitWholeBodyDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input);

    // Destructor
    ~DigitWholeBodyDynamicsConstraints() = default;

    // class methods:
        // return the index of id th dependent joint
    virtual int return_dependent_joint_index(const int id) final override;

        // return the index of id th independent joint
    virtual int return_independent_joint_index(const int id) final override;

        // fill in dependent indeces in a vector
    virtual void fill_dependent_vector(VecX& r, const VecX& v, const bool setZero = false) final override;

        // fill in independent indeces in a vector
    virtual void fill_independent_vector(VecX& r, const VecX& v, const bool setZero = false) final override;

        // fill in dependent columns in a matrix
    virtual void fill_dependent_columns(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in independent columns in a matrix
    virtual void fill_independent_columns(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in dependent rows in a matrix
    virtual void fill_dependent_rows(MatX& r, const MatX& m, const bool setZero = false) final override;

        // fill in independent rows in a matrix
    virtual void fill_independent_rows(MatX& r, const MatX& m, const bool setZero = false) final override;

        // return dependent indeces in a vector
    virtual VecX get_dependent_vector(const VecX& v) final override;

        // return independent indeces in a vector
    virtual VecX get_independent_vector(const VecX& v) final override;

        // return dependent columns in a matrix
    virtual void get_dependent_columns(MatX& r, const MatX& m) final override;

        // return independent columns in a matrix
    virtual void get_independent_columns(MatX& r, const MatX& m) final override;

        // return dependent rows in a matrix
    virtual void get_dependent_rows(MatX& r, const MatX& m) final override;

        // return independent rows in a matrix
    virtual void get_independent_rows(MatX& r, const MatX& m) final override;

        // fill in dependent joint positions in the full joint vector q
        // that satisfies the constraints
        // This usually involves solving inverse kinematics. 
        // You need to implement this method in your derived class!!!
    virtual void setupJointPosition(VecX& q, bool compute_derivatives = true) override {};

        // constraint c(q)
    virtual void get_c(const VecX& q) final override {};

        // J = jacobian(c, q)
    virtual void get_J(const VecX& q) final override;
    
        // Jx_partial_dq = jacobian(J * x, q)
        // where x is a vector that is not dependent on q
    virtual void get_Jx_partial_dq(const VecX& q, const VecX& x) final override {};

        // JTx_partial_dq = jacobian(J^T * x, q)
        // where x is a vector that is not dependent on q
    virtual void get_JTx_partial_dq(const VecX& q, const VecX& x) final override {};

        // Jxy_partial_dq = jacobian(jacobian(J * x, q) * y, q)
        // where x and y are vectors that are not dependent on q
    virtual void get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) final override {};

    // class members:
    std::shared_ptr<Model> modelPtr_ = nullptr;
    std::shared_ptr<Data> dataPtr_ = nullptr;

        // dep/indep joint indeces
    int dependentJointIds[NUM_DEPENDENT_JOINTS] = {0};
    int independentJointIds[NUM_INDEPENDENT_JOINTS] = {0};

        // for contact constraints (forward kinematics mainly)
    pinocchio::FrameIndex left_foot_idx;
    pinocchio::FrameIndex right_foot_idx;

    //     // closed loop related transforms
    //     // you can find all these numbers in digit-v3.xml
    // const Transform left_toeA_rod_endT = Transform(Vec3(0.17 * 2, 0, 0));
    // const Transform left_toeA_anchor_endT = Transform(Vec3(0.0179, -0.009551, -0.054164));

    // const Transform left_toeB_rod_endT = Transform(Vec3(0.144 * 2, 0, 0));
    // const Transform left_toeB_anchor_endT = Transform(Vec3(-0.0181, -0.009551, -0.054164));

    // const Transform left_knee_rod_endT = Transform(Vec3(0.25 * 2, 0, 0));
    // const Transform left_knee_anchor_endT = Transform(Utils::deg2rad(Vec3(4.47, 0.32, 155.8)), 
    //                                                   Vec3(-0.01766, -0.029456, 0.00104)) * 
    //                                         Transform(Vec3(0.113789, -0.011056, 0));

    // const Transform right_toeA_rod_endT = Transform(Vec3(0.17 * 2, 0, 0));
    // const Transform right_toeA_anchor_endT = Transform(Vec3(0.0179, 0.009551, -0.054164));

    // const Transform right_toeB_rod_endT = Transform(Vec3(0.144 * 2, 0, 0));
    // const Transform right_toeB_anchor_endT = Transform(Vec3(-0.0181, 0.009551, -0.054164));
    
    // const Transform right_knee_rod_endT = Transform(Vec3(0.25 * 2, 0, 0));
    // const Transform right_knee_anchor_endT = Transform(Utils::deg2rad(Vec3(-4.47, 0.32, -155.8)), 
    //                                                    Vec3(-0.01766, 0.029456, 0.00104)) *
    //                                          Transform(Vec3(0.113789, 0.011056, 0));

    MatX J_left_foot;
    MatX J_right_foot;
};

}; // namespace DigitWholeBodySysID
}; // namespace RAPTOR

#endif // DIGIT_WHOLEBODY_CONSTRAINTS_H
