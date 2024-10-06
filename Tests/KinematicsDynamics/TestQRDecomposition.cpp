#include "QRDecompositionSolver.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace RAPTOR;

int main() {
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero(); 

    QRDecompositionSolver qrSolver(model); 
    qrSolver.generateRandomObservation();
    qrSolver.computeRegroupMatrix();

    Eigen::VectorXd phi1(qrSolver.phi.size());
    phi1 << qrSolver.beta, qrSolver.phi_d;

    std::cout << "ground truth:\n" << qrSolver.phi.transpose() << std::endl;
    std::cout << "recovered:\n" << (qrSolver.Ginv *  phi1).transpose() << std::endl;

    return 0;
}
