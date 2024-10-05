#include "QRDecompositionSolver.h"
#include "SysIDAlgCombineAllCases.h"
#include "KinovaConstants.h"

using namespace RAPTOR;
using namespace Kinova;

int main() {
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;
    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero(); 

    QRDecompositionSolver qrSolver(model); 
    qrSolver.generateRandomObservation();
    qrSolver.computeRegroupMatrix();

    std::cout << "beta:\n" << qrSolver.beta.transpose() << std::endl;

    Eigen::VectorXd phi1(qrSolver.phi.size());
    phi1 << qrSolver.beta, qrSolver.phi_d;

    std::cout << "ground truth:\n" << qrSolver.phi.transpose() << std::endl;
    std::cout << "recovered:\n" << (qrSolver.Ginv *  phi1).transpose() << std::endl;

    return 0;
}
