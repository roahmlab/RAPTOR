#include "PZDynamics.h"

namespace RAPTOR {
namespace Armour {

KinematicsDynamics::KinematicsDynamics(const std::shared_ptr<RobotInfo>& robotInfoPtr_input,
                                       const std::shared_ptr<BezierCurveInterval>& trajPtr_input) :
    robotInfoPtr_(robotInfoPtr_input),
    trajPtr_(trajPtr_input) {
    const pinocchio::ModelTpl<double>& model = robotInfoPtr_->model;

    // pre-allocate memory
    com_arr = PZsparseArray(model.nv, 3);
    mass_nominal_arr = PZsparseArray(model.nv, 1);
    mass_uncertain_arr = PZsparseArray(model.nv, 1);
    I_nominal_arr = PZsparseArray(model.nv, 9);
    I_uncertain_arr = PZsparseArray(model.nv, 9);

    sphere_centers = PZsparseArray(3 * robotInfoPtr_->num_spheres, NUM_TIME_STEPS);
    torque_nom = PZsparseArray(robotInfoPtr_->num_motors, NUM_TIME_STEPS);
    torque_int = PZsparseArray(robotInfoPtr_->num_motors, NUM_TIME_STEPS);

    contact_force_nom = PZsparseArray(3 * (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors), NUM_TIME_STEPS);
    contact_force_int = PZsparseArray(3 * (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors), NUM_TIME_STEPS);

    contact_moment_nom = PZsparseArray(3 * (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors), NUM_TIME_STEPS);
    contact_moment_int = PZsparseArray(3 * (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors), NUM_TIME_STEPS);

    // initialize robot inertial parameters
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint

        const Vec3 com = model.inertias[pinocchio_joint_id].lever();
        const double mass = model.inertias[pinocchio_joint_id].mass();
        const Mat3 inertia = model.inertias[pinocchio_joint_id].inertia().matrix();

        com_arr(i, 0) = PZsparse(com(0));
        com_arr(i, 1) = PZsparse(com(1));
        com_arr(i, 2) = PZsparse(com(2));

        mass_nominal_arr(i) = PZsparse(mass);
        mass_uncertain_arr(i) = PZsparse(mass, robotInfoPtr_->mass_uncertainty[i]);

        for (int j = 0; j < 3; j++) {
            for (int p = 0; p < 3; p++) {
                I_nominal_arr(i, j * 3 + p) = PZsparse(inertia(p, j));
                I_uncertain_arr(i, j * 3 + p) = PZsparse(inertia(p, j), robotInfoPtr_->inertia_uncertainty[i]);
            }
        }
    }

    // initialize sphere radii
    for (int i = 0; i < robotInfoPtr_->num_spheres; i++) {
        const std::string link_name = "link" + std::to_string(i + 1);

        if (robotInfoPtr_->collision_spheres.find(link_name) == robotInfoPtr_->collision_spheres.end()) {
            continue;
        }

        const auto& spheres = robotInfoPtr_->collision_spheres.at(link_name);
        for (const auto it : spheres) {
            sphere_radii.push_back(it.second);
        }
    }
}

void KinematicsDynamics::reset_trajectory(const std::shared_ptr<BezierCurveInterval>& trajPtr_input) {
    trajPtr_ = trajPtr_input;
}

void KinematicsDynamics::fk(const size_t s_ind) {
    const pinocchio::ModelTpl<double>& model = robotInfoPtr_->model;
    const auto& collision_spheres = robotInfoPtr_->collision_spheres;
    
    // The rotational part of the transformation matrix
    // [[R1, R2, R3],
    //  [R4, R5, R6],
    //  [R7, R8, R9]]
    PZsparse R1(1);
    PZsparse R2(0);
    PZsparse R3(0);
    PZsparse R4(0);
    PZsparse R5(1);
    PZsparse R6(0);
    PZsparse R7(0);
    PZsparse R8(0);
    PZsparse R9(1);

    PZsparse R1_new(1);
    PZsparse R2_new(0);
    PZsparse R3_new(0);
    PZsparse R4_new(0);
    PZsparse R5_new(1);
    PZsparse R6_new(0);
    PZsparse R7_new(0);
    PZsparse R8_new(0);
    PZsparse R9_new(1);

    // The translational part of the transformation matrix
    PZsparse T1(0);
    PZsparse T2(0);
    PZsparse T3(0);

    PZsparse T1_new(0);
    PZsparse T2_new(0);
    PZsparse T3_new(0);

    PZsparseArray R(model.nv, 9); // stored in the following order:
                                  // [[0, 1, 2]
                                  //  [3, 4, 5]
                                  //  [6, 7, 8]]

    // initialize rotation matrices
    for (int i = 0; i < model.nv; i++) {
        const Mat3& jointRotation = model.jointPlacements[i + 1].rotation();

        if (i >= NUM_FACTORS) { // fixed joints, no other rotation
            R(i, 0) = jointRotation(0, 0);
            R(i, 1) = jointRotation(0, 1);
            R(i, 2) = jointRotation(0, 2);
            R(i, 3) = jointRotation(1, 0);
            R(i, 4) = jointRotation(1, 1);
            R(i, 5) = jointRotation(1, 2);
            R(i, 6) = jointRotation(2, 0);
            R(i, 7) = jointRotation(2, 1);
            R(i, 8) = jointRotation(2, 2);

            continue;
        }

        const PZsparse& cq = trajPtr_->cos_q_des(i, s_ind);
        const PZsparse& sq = trajPtr_->sin_q_des(i, s_ind);
        const std::string jointName = model.joints[i + 1].shortname();
        
        if (jointName.find('R') != std::string::npos) { // revolute joint
            if (jointName.find('X') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0);
                R(i, 1) = jointRotation(0, 1) * cq + jointRotation(0, 2) * sq;
                R(i, 2) = -jointRotation(0, 1) * sq + jointRotation(0, 2) * cq;
                R(i, 3) = jointRotation(1, 0);
                R(i, 4) = jointRotation(1, 1) * cq + jointRotation(1, 2) * sq;
                R(i, 5) = -jointRotation(1, 1) * sq + jointRotation(1, 2) * cq;
                R(i, 6) = jointRotation(2, 0);
                R(i, 7) = jointRotation(2, 1) * cq + jointRotation(2, 2) * sq;
                R(i, 8) = -jointRotation(2, 1) * sq + jointRotation(2, 2) * cq;
            }
            else if (jointName.find('Y') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0) * cq - jointRotation(0, 2) * sq;
                R(i, 1) = jointRotation(0, 1);
                R(i, 2) = jointRotation(0, 0) * sq + jointRotation(0, 2) * cq;
                R(i, 3) = jointRotation(1, 0) * cq - jointRotation(1, 2) * sq;
                R(i, 4) = jointRotation(1, 1);
                R(i, 5) = jointRotation(1, 0) * sq + jointRotation(1, 2) * cq;
                R(i, 6) = jointRotation(2, 0) * cq - jointRotation(2, 2) * sq;
                R(i, 7) = jointRotation(2, 1);
                R(i, 8) = jointRotation(2, 0) * sq + jointRotation(2, 2) * cq;
            }
            else if (jointName.find('Z') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0) * cq + jointRotation(0, 1) * sq;
                R(i, 1) = -jointRotation(0, 0) * sq + jointRotation(0, 1) * cq;
                R(i, 2) = jointRotation(0, 2);
                R(i, 3) = jointRotation(1, 0) * cq + jointRotation(1, 1) * sq;
                R(i, 4) = -jointRotation(1, 0) * sq + jointRotation(1, 1) * cq;
                R(i, 5) = jointRotation(1, 2);
                R(i, 6) = jointRotation(2, 0) * cq + jointRotation(2, 1) * sq;
                R(i, 7) = -jointRotation(2, 0) * sq + jointRotation(2, 1) * cq;
                R(i, 8) = jointRotation(2, 2);
            }
            else {
                throw std::invalid_argument("invalid joint type!");
            }
        }
        else if (jointName.find('P') != std::string::npos) { // prismatic joint
            // if (jointName.find('X') != std::string::npos) {
            //     
            // }
            // else if (jointName.find('Y') != std::string::npos) {
            //     
            // }
            // else if (jointName.find('Z') != std::string::npos) {
            //     
            // }
            // else {
                throw std::invalid_argument("Do not support prismatic joints yet!");
            // }
        }
        else {
            throw std::invalid_argument("invalid joint type!");
        }
    }

    // forward kinematics
    int sphere_id = 0;
    for (int i = 0; i < model.nv; i++) {
        const Vec3 jointTranslation = model.jointPlacements[i + 1].translation();

        R1_new = R1 * R(i, 0) + R2 * R(i, 3) + R3 * R(i, 6);
        R2_new = R1 * R(i, 1) + R2 * R(i, 4) + R3 * R(i, 7);
        R3_new = R1 * R(i, 2) + R2 * R(i, 5) + R3 * R(i, 8);
        R4_new = R4 * R(i, 0) + R5 * R(i, 3) + R6 * R(i, 6);
        R5_new = R4 * R(i, 1) + R5 * R(i, 4) + R6 * R(i, 7);
        R6_new = R4 * R(i, 2) + R5 * R(i, 5) + R6 * R(i, 8);
        R7_new = R7 * R(i, 0) + R8 * R(i, 3) + R9 * R(i, 6);
        R8_new = R7 * R(i, 1) + R8 * R(i, 4) + R9 * R(i, 7);
        R9_new = R7 * R(i, 2) + R8 * R(i, 5) + R9 * R(i, 8);

        T1_new = R1 * jointTranslation(0) + R2 * jointTranslation(1) + R3 * jointTranslation(2) + T1;
        T2_new = R4 * jointTranslation(0) + R5 * jointTranslation(1) + R6 * jointTranslation(2) + T2;
        T3_new = R7 * jointTranslation(0) + R8 * jointTranslation(1) + R9 * jointTranslation(2) + T3;

        R1 = R1_new;
        R2 = R2_new;
        R3 = R3_new;
        R4 = R4_new;
        R5 = R5_new;
        R6 = R6_new;
        R7 = R7_new;
        R8 = R8_new;
        R9 = R9_new;

        T1 = T1_new;
        T2 = T2_new;
        T3 = T3_new;

        std::string current_link_name = "link" + std::to_string(i + 1);
        if (collision_spheres.find(current_link_name) != collision_spheres.end()) {
            const auto& spheres = collision_spheres.at(current_link_name);
            for (const auto it : spheres) {
                const Vec3 sphere_offset = it.first;

                T1_new = R1 * sphere_offset(0) + R2 * sphere_offset(1) + R3 * sphere_offset(2) + T1;
                T2_new = R4 * sphere_offset(0) + R5 * sphere_offset(1) + R6 * sphere_offset(2) + T2;
                T3_new = R7 * sphere_offset(0) + R8 * sphere_offset(1) + R9 * sphere_offset(2) + T3;

                sphere_centers(sphere_id * 3 + 0, s_ind) = T1_new;
                sphere_centers(sphere_id * 3 + 1, s_ind) = T2_new;
                sphere_centers(sphere_id * 3 + 2, s_ind) = T3_new;
                sphere_id++;
            }
        }
    }
}

void KinematicsDynamics::rnea(const size_t s_ind,
                              const PZsparseArray& mass_arr,
                              const PZsparseArray& I_arr,
                              PZsparseArray& u,
                              PZsparseArray& contact_force,
			                  PZsparseArray& contact_moment) {
    const pinocchio::ModelTpl<double>& model = robotInfoPtr_->model;

    // the following terms are automatically initialized to zero
    PZsparse w1;
    PZsparse w2;
    PZsparse w3;
    PZsparse wdot1;
    PZsparse wdot2;
    PZsparse wdot3;
    PZsparse w_aux1;
    PZsparse w_aux2;
    PZsparse w_aux3;
    PZsparse linear_acc1 = -model.gravity.linear()(0);
    PZsparse linear_acc2 = -model.gravity.linear()(1);
    PZsparse linear_acc3 = -model.gravity.linear()(2);

    PZsparse temp1;
    PZsparse temp2;
    PZsparse temp3;
    PZsparse temp4;
    PZsparse temp5;
    PZsparse temp6;

    PZsparseArray F(model.nv, 3);
    PZsparseArray N(model.nv, 3);
    PZsparseArray R(model.nv, 9); // stored in the following order:
                                  // [[0, 1, 2]
                                  //  [3, 4, 5]
                                  //  [6, 7, 8]]

    // initialize rotation matrices
    for (int i = 0; i < model.nv; i++) {
        const Mat3& jointRotation = model.jointPlacements[i + 1].rotation();
        
        if (i >= NUM_FACTORS) { // fixed joints, no other rotation
            R(i, 0) = jointRotation(0, 0);
            R(i, 1) = jointRotation(0, 1);
            R(i, 2) = jointRotation(0, 2);
            R(i, 3) = jointRotation(1, 0);
            R(i, 4) = jointRotation(1, 1);
            R(i, 5) = jointRotation(1, 2);
            R(i, 6) = jointRotation(2, 0);
            R(i, 7) = jointRotation(2, 1);
            R(i, 8) = jointRotation(2, 2);

            continue;
        }

        const PZsparse& cq = trajPtr_->cos_q_des(i, s_ind);
        const PZsparse& sq = trajPtr_->sin_q_des(i, s_ind);
        const std::string jointName = model.joints[i + 1].shortname();

        if (jointName.find('R') != std::string::npos) { // revolute joint
            if (jointName.find('X') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0);
                R(i, 1) = jointRotation(0, 1) * cq + jointRotation(0, 2) * sq;
                R(i, 2) = -jointRotation(0, 1) * sq + jointRotation(0, 2) * cq;
                R(i, 3) = jointRotation(1, 0);
                R(i, 4) = jointRotation(1, 1) * cq + jointRotation(1, 2) * sq;
                R(i, 5) = -jointRotation(1, 1) * sq + jointRotation(1, 2) * cq;
                R(i, 6) = jointRotation(2, 0);
                R(i, 7) = jointRotation(2, 1) * cq + jointRotation(2, 2) * sq;
                R(i, 8) = -jointRotation(2, 1) * sq + jointRotation(2, 2) * cq;
            }
            else if (jointName.find('Y') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0) * cq - jointRotation(0, 2) * sq;
                R(i, 1) = jointRotation(0, 1);
                R(i, 2) = jointRotation(0, 0) * sq + jointRotation(0, 2) * cq;
                R(i, 3) = jointRotation(1, 0) * cq - jointRotation(1, 2) * sq;
                R(i, 4) = jointRotation(1, 1);
                R(i, 5) = jointRotation(1, 0) * sq + jointRotation(1, 2) * cq;
                R(i, 6) = jointRotation(2, 0) * cq - jointRotation(2, 2) * sq;
                R(i, 7) = jointRotation(2, 1);
                R(i, 8) = jointRotation(2, 0) * sq + jointRotation(2, 2) * cq;
            }
            else if (jointName.find('Z') != std::string::npos) {
                R(i, 0) = jointRotation(0, 0) * cq + jointRotation(0, 1) * sq;
                R(i, 1) = -jointRotation(0, 0) * sq + jointRotation(0, 1) * cq;
                R(i, 2) = jointRotation(0, 2);
                R(i, 3) = jointRotation(1, 0) * cq + jointRotation(1, 1) * sq;
                R(i, 4) = -jointRotation(1, 0) * sq + jointRotation(1, 1) * cq;
                R(i, 5) = jointRotation(1, 2);
                R(i, 6) = jointRotation(2, 0) * cq + jointRotation(2, 1) * sq;
                R(i, 7) = -jointRotation(2, 0) * sq + jointRotation(2, 1) * cq;
                R(i, 8) = jointRotation(2, 2);
            }
            else {
                throw std::invalid_argument("invalid joint type!");
            }
        }
        else if (jointName.find('P') != std::string::npos) { // prismatic joint
            // if (jointName.find('X') != std::string::npos) {
            //     
            // }
            // else if (jointName.find('Y') != std::string::npos) {
            //     
            // }
            // else if (jointName.find('Z') != std::string::npos) {
            //     
            // }
            // else {
                throw std::invalid_argument("Do not support prismatic joints yet!");
            // }
        }
        else {
            throw std::invalid_argument("invalid joint type!");
        }   
    }

    // RNEA forward recursion
    for (int i = 0; i < model.nv; i++) {
        const PZsparse& cq = trajPtr_->cos_q_des(i, s_ind);
        const PZsparse& sq = trajPtr_->sin_q_des(i, s_ind);
        const PZsparse& qd = trajPtr_->qd_des(i, s_ind);
        const PZsparse& qda = trajPtr_->qda_des(i, s_ind);
        const PZsparse& qdda = trajPtr_->qdda_des(i, s_ind);
        const Vec3 jointTranslation = model.jointPlacements[i + 1].translation();

        // NOTE:
        // This is just a simplified implementation!!!
        // We assume all fixed joints are at the end and the revolute joints are consecutive
        if (i < robotInfoPtr_->num_motors) { // revolute joints
            // line 16
            // linear_acc = trajPtr_->R_t(i, s_ind) * (linear_acc 
            //                                      + cross(wdot, trans_matrix(i, 0)) 
            //                                      + cross(w, cross(w_aux, trans_matrix(i, 0))));
            temp1 = w_aux2 * jointTranslation(2) - w_aux3 * jointTranslation(1);
            temp2 = w_aux3 * jointTranslation(0) - w_aux1 * jointTranslation(2);
            temp3 = w_aux1 * jointTranslation(1) - w_aux2 * jointTranslation(0);
            temp4 = w2 * temp3 - w3 * temp2;
            temp5 = w3 * temp1 - w1 * temp3;
            temp6 = w1 * temp2 - w2 * temp1;
            temp1 = wdot2 * jointTranslation(2) - wdot3 * jointTranslation(1);
            temp2 = wdot3 * jointTranslation(0) - wdot1 * jointTranslation(2);
            temp3 = wdot1 * jointTranslation(1) - wdot2 * jointTranslation(0);
            temp1 = linear_acc1 + temp1 + temp4;
            temp2 = linear_acc2 + temp2 + temp5;
            temp3 = linear_acc3 + temp3 + temp6;
            linear_acc1 = R(i, 0) * temp1 + R(i, 3) * temp2 + R(i, 6) * temp3;
            linear_acc2 = R(i, 1) * temp1 + R(i, 4) * temp2 + R(i, 7) * temp3;
            linear_acc3 = R(i, 2) * temp1 + R(i, 5) * temp2 + R(i, 8) * temp3;

            // line 13
            // w = trajPtr_->R_t(i, s_ind) * w;
            // w.addOneDimPZ(trajPtr_->qd_des(i, s_ind), abs(axes[i]) - 1, 0);
            temp1 = R(i, 0) * w1 + R(i, 3) * w2 + R(i, 6) * w3;
            temp2 = R(i, 1) * w1 + R(i, 4) * w2 + R(i, 7) * w3;
            temp3 = R(i, 2) * w1 + R(i, 5) * w2 + R(i, 8) * w3;
            w1 = temp1;
            w2 = temp2;
            w3 = temp3;
            applyOnlyOneDimension(model.joints[i + 1].shortname(), 
                                  qd,
                                  w1, w2, w3);

            // line 14
            // w_aux = trajPtr_->R_t(i, s_ind) * w_aux;
            temp1 = R(i, 0) * w_aux1 + R(i, 3) * w_aux2 + R(i, 6) * w_aux3;
            temp2 = R(i, 1) * w_aux1 + R(i, 4) * w_aux2 + R(i, 7) * w_aux3;
            temp3 = R(i, 2) * w_aux1 + R(i, 5) * w_aux2 + R(i, 8) * w_aux3;
            w_aux1 = temp1;
            w_aux2 = temp2;
            w_aux3 = temp3;

            // line 15
            // wdot = trajPtr_->R_t(i, s_ind) * wdot;
            temp1 = R(i, 0) * wdot1 + R(i, 3) * wdot2 + R(i, 6) * wdot3;
            temp2 = R(i, 1) * wdot1 + R(i, 4) * wdot2 + R(i, 7) * wdot3;
            temp3 = R(i, 2) * wdot1 + R(i, 5) * wdot2 + R(i, 8) * wdot3;
            wdot1 = temp1;
            wdot2 = temp2;
            wdot3 = temp3;

            // PZsparse temp(3, 1); // temp = joint_vel(robot_params.q_index(i))*z(:,i)
            // temp.addOneDimPZ(trajPtr_->qd_des(i, s_ind), abs(axes[i]) - 1, 0);
            // wdot = wdot + cross(w_aux, temp);
            crossOnlyOneDimension(model.joints[i + 1].shortname(), 
                                  qd,
                                  w_aux1, w_aux2, w_aux3,
                                  temp1, temp2, temp3);
            wdot1 += temp1;
            wdot2 += temp2;
            wdot3 += temp3;

            // wdot.addOneDimPZ(trajPtr_->qdda_des(i, s_ind), abs(axes[i]) - 1, 0);
            applyOnlyOneDimension(model.joints[i + 1].shortname(), 
                                  qdda,
                                  wdot1, wdot2, wdot3);

            // line 14
            // w_aux.addOneDimPZ(trajPtr_->qda_des(i, s_ind), abs(axes[i]) - 1, 0);
            applyOnlyOneDimension(model.joints[i + 1].shortname(), 
                                  qda,
                                  w_aux1, w_aux2, w_aux3);
        }
        else { // fixed joints
            // line 16
            // linear_acc = trajPtr_->R_t(i, s_ind) * (linear_acc 
            //                                      + cross(wdot, trans_matrix(i, 0)) 
            //                                      + cross(w, cross(w_aux, trans_matrix(i, 0))));
            temp1 = w_aux2 * jointTranslation(2) - w_aux3 * jointTranslation(1);
            temp2 = w_aux3 * jointTranslation(0) - w_aux1 * jointTranslation(2);
            temp3 = w_aux1 * jointTranslation(1) - w_aux2 * jointTranslation(0);
            temp4 = w2 * temp3 - w3 * temp2;
            temp5 = w3 * temp1 - w1 * temp3;
            temp6 = w1 * temp2 - w2 * temp1;
            temp1 = wdot2 * jointTranslation(2) - wdot3 * jointTranslation(1);
            temp2 = wdot3 * jointTranslation(0) - wdot1 * jointTranslation(2);
            temp3 = wdot1 * jointTranslation(1) - wdot2 * jointTranslation(0);
            temp1 = linear_acc1 + temp1 + temp4;
            temp2 = linear_acc2 + temp2 + temp5;
            temp3 = linear_acc3 + temp3 + temp6;
            linear_acc1 = R(i, 0) * temp1 + R(i, 3) * temp2 + R(i, 6) * temp3;
            linear_acc2 = R(i, 1) * temp1 + R(i, 4) * temp2 + R(i, 7) * temp3;
            linear_acc3 = R(i, 2) * temp1 + R(i, 5) * temp2 + R(i, 8) * temp3;

            // line 13
            // w = trajPtr_->R_t(i, s_ind) * w;
            temp1 = R(i, 0) * temp1 + R(i, 3) * temp2 + R(i, 6) * temp3;
            temp2 = R(i, 1) * temp1 + R(i, 4) * temp2 + R(i, 7) * temp3;
            temp3 = R(i, 2) * temp1 + R(i, 5) * temp2 + R(i, 8) * temp3;
            w1 = temp1;
            w2 = temp2;
            w3 = temp3;

            // line 14
            // w_aux = trajPtr_->R_t(i, s_ind) * w_aux;
            temp1 = R(i, 0) * w_aux1 + R(i, 3) * w_aux2 + R(i, 6) * w_aux3;
            temp2 = R(i, 1) * w_aux1 + R(i, 4) * w_aux2 + R(i, 7) * w_aux3;
            temp3 = R(i, 2) * w_aux1 + R(i, 5) * w_aux2 + R(i, 8) * w_aux3;
            w_aux1 = temp1;
            w_aux2 = temp2;
            w_aux3 = temp3;

            // line 15
            // wdot = trajPtr_->R_t(i, s_ind) * wdot;
            temp1 = R(i, 0) * wdot1 + R(i, 3) * wdot2 + R(i, 6) * wdot3;
            temp2 = R(i, 1) * wdot1 + R(i, 4) * wdot2 + R(i, 7) * wdot3;
            temp3 = R(i, 2) * wdot1 + R(i, 5) * wdot2 + R(i, 8) * wdot3;
            wdot1 = temp1;
            wdot2 = temp2;
            wdot3 = temp3;
        }

        // line 23 & 27
        // F(i, 0) = mass_arr(i, 0) * (linear_acc
        //                              + cross(wdot, com_matrix(i, 0))
        //                              + cross(w, cross(w_aux, com_matrix(i, 0))));
        temp1 = w_aux2 * com_arr(i, 2) - w_aux3 * com_arr(i, 1);
        temp2 = w_aux3 * com_arr(i, 0) - w_aux1 * com_arr(i, 2);
        temp3 = w_aux1 * com_arr(i, 1) - w_aux2 * com_arr(i, 0);
        temp4 = w2 * temp3 - w3 * temp2;
        temp5 = w3 * temp1 - w1 * temp3;
        temp6 = w1 * temp2 - w2 * temp1;
        temp1 = wdot2 * com_arr(i, 2) - wdot3 * com_arr(i, 1);
        temp2 = wdot3 * com_arr(i, 0) - wdot1 * com_arr(i, 2);
        temp3 = wdot1 * com_arr(i, 1) - wdot2 * com_arr(i, 0);
        F(i, 0) = mass_arr(i) * (linear_acc1 + temp1 + temp4);
        F(i, 1) = mass_arr(i) * (linear_acc2 + temp2 + temp5);
        F(i, 2) = mass_arr(i) * (linear_acc3 + temp3 + temp6);

        // line 29
        // N(i, 0) = I_arr(i, 0) * wdot + cross(w_aux, (I_arr(i, 0) * w));
        temp1 = I_arr(i, 0) * w1 + I_arr(i, 1) * w2 + I_arr(i, 2) * w3;
        temp2 = I_arr(i, 3) * w1 + I_arr(i, 4) * w2 + I_arr(i, 5) * w3;
        temp3 = I_arr(i, 6) * w1 + I_arr(i, 7) * w2 + I_arr(i, 8) * w3;
        temp4 = w_aux2 * temp3 - w_aux3 * temp2;
        temp5 = w_aux3 * temp1 - w_aux1 * temp3;
        temp6 = w_aux1 * temp2 - w_aux2 * temp1;
        temp1 = I_arr(i, 0) * wdot1 + I_arr(i, 1) * wdot2 + I_arr(i, 2) * wdot3;
        temp2 = I_arr(i, 3) * wdot1 + I_arr(i, 4) * wdot2 + I_arr(i, 5) * wdot3;
        temp3 = I_arr(i, 6) * wdot1 + I_arr(i, 7) * wdot2 + I_arr(i, 8) * wdot3;
        N(i, 0) = temp1 + temp4;
        N(i, 1) = temp2 + temp5;
        N(i, 2) = temp3 + temp6;
    }

    // RNEA reverse recursion
    for (int i = model.nv - 1; i >= 0; i--) {
        const PZsparse& cq = trajPtr_->cos_q_des(i, s_ind);
        const PZsparse& sq = trajPtr_->sin_q_des(i, s_ind);
        const PZsparse& qd = trajPtr_->qd_des(i, s_ind);
        const PZsparse& qda = trajPtr_->qda_des(i, s_ind);
        const PZsparse& qdda = trajPtr_->qdda_des(i, s_ind);
        const std::string jointName = model.joints[i + 1].shortname();

        if (i == model.nv - 1) {
            // line 29
            // n = N(i, 0)
            //     + cross(com_matrix(i, 0), F(i, 0));
            temp1 = com_arr(i, 1) * F(i, 2) - com_arr(i, 2) * F(i, 1);
            temp2 = com_arr(i, 2) * F(i, 0) - com_arr(i, 0) * F(i, 2);
            temp3 = com_arr(i, 0) * F(i, 1) - com_arr(i, 1) * F(i, 0);
            N(i, 0) += temp1;
            N(i, 1) += temp2;
            N(i, 2) += temp3;

            // line 28
            // f = F(i, 0);
        }
        else {
            const Vec3 parentJointTranslation = model.jointPlacements[i + 1 + 1].translation();
            
            // line 29
            // n = N(i, 0)
            //     + trajPtr_->R(i + 1, s_ind) * n
            //     + cross(com_matrix(i, 0), F(i, 0))
            //     + cross(trans_matrix(i + 1, 0), trajPtr_->R(i + 1, s_ind) * f);
            temp1 = R(i + 1, 0) * F(i + 1, 0) + R(i + 1, 1) * F(i + 1, 1) + R(i + 1, 2) * F(i + 1, 2);
            temp2 = R(i + 1, 3) * F(i + 1, 0) + R(i + 1, 4) * F(i + 1, 1) + R(i + 1, 5) * F(i + 1, 2);
            temp3 = R(i + 1, 6) * F(i + 1, 0) + R(i + 1, 7) * F(i + 1, 1) + R(i + 1, 8) * F(i + 1, 2);
            temp4 = parentJointTranslation(1) * temp3 - parentJointTranslation(2) * temp2;
            temp5 = parentJointTranslation(2) * temp1 - parentJointTranslation(0) * temp3;
            temp6 = parentJointTranslation(0) * temp2 - parentJointTranslation(1) * temp1;
            temp1 = temp4 + com_arr(i, 1) * F(i, 2) - com_arr(i, 2) * F(i, 1);
            temp2 = temp5 + com_arr(i, 2) * F(i, 0) - com_arr(i, 0) * F(i, 2);
            temp3 = temp6 + com_arr(i, 0) * F(i, 1) - com_arr(i, 1) * F(i, 0);
            temp4 = R(i + 1, 0) * N(i + 1, 0) + R(i + 1, 1) * N(i + 1, 1) + R(i + 1, 2) * N(i + 1, 2);
            temp5 = R(i + 1, 3) * N(i + 1, 0) + R(i + 1, 4) * N(i + 1, 1) + R(i + 1, 5) * N(i + 1, 2);
            temp6 = R(i + 1, 6) * N(i + 1, 0) + R(i + 1, 7) * N(i + 1, 1) + R(i + 1, 8) * N(i + 1, 2);
            N(i, 0) += temp1 + temp4;
            N(i, 1) += temp2 + temp5;
            N(i, 2) += temp3 + temp6;

            // line 28
            // f = trajPtr_->R(i + 1, s_ind) * f + F(i, 0);
            temp1 = R(i + 1, 0) * F(i + 1, 0) + R(i + 1, 1) * F(i + 1, 1) + R(i + 1, 2) * F(i + 1, 2);
            temp2 = R(i + 1, 3) * F(i + 1, 0) + R(i + 1, 4) * F(i + 1, 1) + R(i + 1, 5) * F(i + 1, 2);
            temp3 = R(i + 1, 6) * F(i + 1, 0) + R(i + 1, 7) * F(i + 1, 1) + R(i + 1, 8) * F(i + 1, 2);
            F(i, 0) += temp1;
            F(i, 1) += temp2;
            F(i, 2) += temp3;
        }

        if (i < robotInfoPtr_->num_motors) {
            if (jointName.find('R') != std::string::npos) { // revolute joint
                if (jointName.find('X') != std::string::npos) {
                    u(i, s_ind) = N(i, 0);
                }
                else if (jointName.find('Y') != std::string::npos) {
                    u(i, s_ind) = N(i, 1);
                }
                else if (jointName.find('Z') != std::string::npos) {
                    u(i, s_ind) = N(i, 2);
                }
                else {
                    throw std::invalid_argument("robotInfoPtr_->num_motorsinvalid joint type!");
                }
            }
            else if (jointName.find('P') != std::string::npos) { // prismatic joint
                // if (jointName.find('X') != std::string::npos) {
                //     
                // }
                // else if (jointName.find('Y') != std::string::npos) {
                //     
                // }
                // else if (jointName.find('Z') != std::string::npos) {
                //     
                // }
                // else {
                    throw std::invalid_argument("Do not support prismatic joints yet!");
                // }
            }
            else {
                throw std::invalid_argument("robotInfoPtr_->num_motorsinvalid joint type!");
            }

            u(i, s_ind) += model.damping[i] * qd + model.armature[i] * qdda;

            // friction is directly applied on the torque limits
        }
        else {
            contact_force(3 * (i - robotInfoPtr_->num_motors) + 0, s_ind) = F(i, 0);
            contact_force(3 * (i - robotInfoPtr_->num_motors) + 1, s_ind) = F(i, 1);
            contact_force(3 * (i - robotInfoPtr_->num_motors) + 2, s_ind) = F(i, 2);

            contact_moment(3 * (i - robotInfoPtr_->num_motors) + 0, s_ind) = N(i, 0);
            contact_moment(3 * (i - robotInfoPtr_->num_motors) + 1, s_ind) = N(i, 1);
            contact_moment(3 * (i - robotInfoPtr_->num_motors) + 2, s_ind) = N(i, 2);
        }
    }
}

void applyOnlyOneDimension(const std::string& jointName,
						   const PZsparse& in,
						   PZsparse& out1,
						   PZsparse& out2,
						   PZsparse& out3) {
    if (jointName.find('R') != std::string::npos) { // revolute joint
        if (jointName.find('X') != std::string::npos) {
            out1 += in;
        }
        else if (jointName.find('Y') != std::string::npos) {
            out2 += in;
        }
        else if (jointName.find('Z') != std::string::npos) {
            out3 += in;
        }
        else {
            throw std::invalid_argument("applyOnlyOneDimension: invalid joint type!");
        }
    }
    else if (jointName.find('P') != std::string::npos) { // prismatic joint
        // if (jointName.find('X') != std::string::npos) {
        //     
        // }
        // else if (jointName.find('Y') != std::string::npos) {
        //     
        // }
        // else if (jointName.find('Z') != std::string::npos) {
        //     
        // }
        // else {
            throw std::invalid_argument("Do not support prismatic joints yet!");
        // }
    }
    else {
        throw std::invalid_argument("applyOnlyOneDimension: invalid joint type!");
    }
}

void crossOnlyOneDimension(const std::string& jointName,
						   const PZsparse& v,
						   const PZsparse& in1,
						   const PZsparse& in2,
						   const PZsparse& in3,
						   PZsparse& out1,
						   PZsparse& out2,
						   PZsparse& out3) {
    if (jointName.find('R') != std::string::npos) { // revolute joint
        if (jointName.find('X') != std::string::npos) {
            out1 = 0;
            out2 = in3 * v;
            out3 = -in2 * v;
        }
        else if (jointName.find('Y') != std::string::npos) {
            out1 = -in3 * v;
            out2 = 0;
            out3 = in1 * v;
        }
        else if (jointName.find('Z') != std::string::npos) {
            out1 = in2 * v;
            out2 = -in1 * v;
            out3 = 0;
        }
        else {
            throw std::invalid_argument("crossOnlyOneDimension: invalid joint type!");
        }
    }
    else if (jointName.find('P') != std::string::npos) { // prismatic joint
        // if (jointName.find('X') != std::string::npos) {
        //     
        // }
        // else if (jointName.find('Y') != std::string::npos) {
        //     
        // }
        // else if (jointName.find('Z') != std::string::npos) {
        //     
        // }
        // else {
            throw std::invalid_argument("Do not support prismatic joints yet!");
        // }
    }
    else {
        throw std::invalid_argument("crossOnlyOneDimension: invalid joint type!");
    }
}

}; // namespace Armour
}; // namespace RAPTOR