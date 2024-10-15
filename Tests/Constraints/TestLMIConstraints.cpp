#include "LMIConstraints.h"

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
    lmi.compute(phi, true, true);

    std::cout << "g:\n" << lmi.g.transpose() << std::endl;

    // // test the gradient
    // const Eigen::MatrixXd& J_analytical = lmi.pg_pz;
    // Eigen::MatrixXd J_numerical = 
    //     Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());

    // for (int i = 0; i < phi.size(); i++) {
    //     Eigen::VectorXd phi_plus = phi;
    //     phi_plus(i) += 1e-7;
    //     lmi.compute(phi_plus, false);
    //     const Eigen::VectorXd g_plus = lmi.g;
    //     Eigen::VectorXd phi_minus = phi;
    //     phi_minus(i) -= 1e-7;
    //     lmi.compute(phi_minus, false);
    //     const Eigen::VectorXd g_minus = lmi.g;
    //     J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    // }
    // std::cout << "jacobian total difference: " << (J_analytical - J_numerical).norm() << std::endl;

    // // test hessian
    // Eigen::Array<Eigen::MatrixXd, 1, Eigen::Dynamic> H_analytical = lmi.pg_pz_pz;
    // for (int i = 0; i < phi.size(); i++) {
    //     Eigen::VectorXd q_plus = phi;
    //     q_plus(i) += 1e-7;
    //     lmi.compute(q_plus, true, false);
    //     const Eigen::MatrixXd J_plus = lmi.pg_pz;
    //     Eigen::VectorXd q_minus = phi;
    //     q_minus(i) -= 1e-7;
    //     lmi.compute(q_minus, true, false);
    //     const Eigen::MatrixXd J_minus = lmi.pg_pz;
    //     const Eigen::MatrixXd H_numerical_row = (J_plus - J_minus) / 2e-7;

    //     for (int j = 0; j < H_analytical.size(); j++) {
    //         std::cout << "hessian row difference: " << (H_analytical(j).row(i) - H_numerical_row.row(j)).norm() << std::endl;
    //     }
    // }

    return 0;
}