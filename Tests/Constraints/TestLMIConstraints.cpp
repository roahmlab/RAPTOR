#include "LMIConstraints.h"
#include "KinematicsConstraints.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iostream>

using namespace RAPTOR;

int main() {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // recover the inertial parameters from the model
    // which is a physically consistent set of parameters
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(10 * model.nv); 
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint
        phi.segment<10>(10 * i) = 
            model.inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }

    LMIConstraints lmi(model.nv, 10 * model.nv);
    lmi.compute(phi, true);

    std::cout << "g:\n" << lmi.g.transpose() << std::endl;

    // test the gradient
    const Eigen::MatrixXd& J_analytical = lmi.pg_pz;
    Eigen::MatrixXd J_numerical = 
        Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());

    for (int i = 0; i < phi.size(); i++) {
        Eigen::VectorXd phi_plus = phi;
        phi_plus(i) += 1e-7;
        lmi.compute(phi_plus, false);
        const Eigen::VectorXd g_plus = lmi.g;
        Eigen::VectorXd phi_minus = phi;
        phi_minus(i) -= 1e-7;
        lmi.compute(phi_minus, false);
        const Eigen::VectorXd g_minus = lmi.g;
        J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    }
    std::cout << "difference: " << (J_analytical - J_numerical).norm() << std::endl;

    return 0;
}