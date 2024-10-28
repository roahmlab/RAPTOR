#include "QRDecompositionSolver.h"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include <pinocchio/algorithm/rnea.hpp>
using namespace RAPTOR;
int main() {
    std::srand(time(nullptr));
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();
    QRDecompositionSolver qrSolver(model);
    qrSolver.generateRandomObservation();
    qrSolver.computeRegroupMatrix();
    Eigen::VectorXd phi1(qrSolver.phi.size());
    phi1 << qrSolver.beta, qrSolver.phi_d;
    // std::cout << "ground truth:\n" << qrSolver.phi.transpose() << std::endl;
    // std::cout << "recovered:\n" << (qrSolver.Ginv * phi1).transpose() << std::endl;
    std::cout << qrSolver.phi.transpose() -
                 (qrSolver.Ginv * phi1).transpose() << std::endl;
    Eigen::VectorXd q = 2* M_PI* Eigen::VectorXd::Random(model.nv);
    Eigen::VectorXd v = 2* M_PI*Eigen::VectorXd::Random(model.nv);
    Eigen::VectorXd a = 2* M_PI*Eigen::VectorXd::Random(model.nv);
    pinocchio::computeJointTorqueRegressor(model, data, q, v, a);
    pinocchio::rnea(model, data, q, v, a);
    std::cout << "rnea tau: "<< (data.tau.transpose()) << std::endl;
    std::cout << "full matrix: "<< (data.jointTorqueRegressor* qrSolver.phi).transpose() << std::endl;
    std::cout << "regroup: "<< (data.jointTorqueRegressor* qrSolver.Aid * qrSolver.phi_id +
            data.jointTorqueRegressor* qrSolver.Ad * qrSolver.phi_d).transpose() << std::endl;
    std::cout << "difference: " << data.jointTorqueRegressor * qrSolver.Aid * qrSolver.Kd - data.jointTorqueRegressor * qrSolver.Ad<< std::endl;
    return 0;
}