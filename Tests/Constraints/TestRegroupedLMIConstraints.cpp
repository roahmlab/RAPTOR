#include "RegroupedLMIConstraints.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iostream>

using namespace RAPTOR;

int main() {
    std::srand(time(nullptr));
    
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero(); 

    // Define QR decomposition solver
    std::shared_ptr<QRDecompositionSolver> qrSolverPtr_ = 
        std::make_shared<QRDecompositionSolver>(model); 
    qrSolverPtr_->generateRandomObservation();
    qrSolverPtr_->computeRegroupMatrix();

    // Define the LMI constraints
    LMIConstraints lmi(model.nv, 10 * model.nv);
    lmi.compute(qrSolverPtr_->phi, false);
    
    // Define the regrouped LMI constraints
    RegroupedLMIConstraints regroupedLMI(qrSolverPtr_, model.nv, 10 * model.nv);
    Eigen::VectorXd z(qrSolverPtr_->dim_id + qrSolverPtr_->dim_d);
    z << qrSolverPtr_->beta, qrSolverPtr_->phi_d;
    regroupedLMI.compute(z, true, true);

    // test the constraints consistency
    std::cout << "difference: " << (lmi.g - regroupedLMI.g).norm() << std::endl;

    // test the gradient
    const Eigen::MatrixXd& J_analytical = regroupedLMI.pg_pz;
    Eigen::MatrixXd J_numerical = 
        Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());

    for (int i = 0; i < z.size(); i++) {
        Eigen::VectorXd z_plus = z;
        z_plus(i) += 1e-7;
        regroupedLMI.compute(z_plus, false);
        const Eigen::VectorXd g_plus = regroupedLMI.g;
        Eigen::VectorXd z_minus = z;
        z_minus(i) -= 1e-7;
        regroupedLMI.compute(z_minus, false);
        const Eigen::VectorXd g_minus = regroupedLMI.g;
        J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    }
    std::cout << "jacobian total difference: " << (J_analytical - J_numerical).norm() << std::endl;

    // test hessian
    Eigen::Array<Eigen::MatrixXd, 1, Eigen::Dynamic> H_analytical = regroupedLMI.pg_pz_pz;
    for (int i = 0; i < z.size(); i++) {
        Eigen::VectorXd q_plus = z;
        q_plus(i) += 1e-7;
        regroupedLMI.compute(q_plus, true, false);
        const Eigen::MatrixXd J_plus = regroupedLMI.pg_pz;
        Eigen::VectorXd q_minus = z;
        q_minus(i) -= 1e-7;
        regroupedLMI.compute(q_minus, true, false);
        const Eigen::MatrixXd J_minus = regroupedLMI.pg_pz;
        const Eigen::MatrixXd H_numerical_row = (J_plus - J_minus) / 2e-7;

        for (int j = 0; j < H_analytical.size(); j++) {
            std::cout << "hessian row difference: " << (H_analytical(j).row(i) - H_numerical_row.row(j)).norm() << std::endl;
        }
    }

    return 0;
}