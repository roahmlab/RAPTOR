#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "RegressorInverseDynamics.h"

using namespace IDTO;
using namespace Kinova;
using namespace Ipopt;

int main() {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);
    
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.81;
    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero();

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    return 0;
}