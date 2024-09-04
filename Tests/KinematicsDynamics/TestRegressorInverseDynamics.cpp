#include "RegressorInverseDynamics.h"
#include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "Trajectories.h"

using VecX = Eigen::VectorXd;
using namespace RAPTOR;

int main() {
    // Define robot model
    
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);
    // Create a trajectory
    int N = 50;  // number of time steps
    double T = 10.0;  // total time
    std::shared_ptr<Trajectories> trajPtr = std::make_shared<Trajectories>(model.nv * 3 * N, T, N, model.nv, Uniform);
    // Initialize RegressorInverseDynamics
    RegressorInverseDynamics regressor_id(model, trajPtr);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    VecX z = M_PI * VecX::Random(trajPtr->varLength);

    // Compute inverse dynamics using RegressorInverseDynamics
    auto start_clock = std::chrono::high_resolution_clock::now();
    regressor_id.compute(z, false);
    auto stop_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "RegressorInverseDynamics: " << duration.count() << " nanoseconds" << std::endl;

    // Compute inverse dynamics using pinocchio::rnea
    start_clock = std::chrono::high_resolution_clock::now();
    VecX tau_pinocchio = pinocchio::rnea(model, data, trajPtr->q(0), trajPtr->q_d(0), trajPtr->q_dd(0));
    stop_clock = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "Pinocchio RNEA: " << duration.count() << " nanoseconds" << std::endl;

    // Compare the results
    std::cout << "RegressorInverseDynamics result:" << std::endl;
    std::cout << regressor_id.tau(0).transpose() << std::endl;
    std::cout << "Pinocchio RNEA result:" << std::endl;
    std::cout << tau_pinocchio.transpose() << std::endl;



    return 0;
}