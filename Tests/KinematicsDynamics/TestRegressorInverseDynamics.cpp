#include "RegressorInverseDynamics.h"
#include <chrono>
#include "pinocchio/algorithm/rnea.hpp"
#include "Polynomials.h"

using VecX = Eigen::VectorXd;
using namespace RAPTOR;

int main() {
    // Define robot model
    
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    pinocchio::ModelTpl<double> model = model_double.cast<double>();
    pinocchio::DataTpl<double> data(model);

    // Disable rotor inertia, friction, and damping
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

    // Create a trajectory
    int N = 5;  // number of time steps
    double T = 10.0;  // total time
    int degree = 5;  // degree of the polynomial
    std::shared_ptr<Trajectories> trajPtr = 
        std::make_shared<Polynomials>(T, N, model.nv, Uniform, degree);

    // Initialize RegressorInverseDynamics
    RegressorInverseDynamics regressor_id(model, trajPtr);

    // Generate random joint p, v, and a (not accurate)
    std::srand(std::time(nullptr));
    VecX z = M_2_PI * VecX::Random(trajPtr->varLength).array() - M_PI;

    // Compute inverse dynamics using RegressorInverseDynamics
    auto start_clock = std::chrono::high_resolution_clock::now();
    regressor_id.compute(z, false);
    auto stop_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "RegressorInverseDynamics: " << duration.count() << " nanoseconds" << std::endl;

    // Compute inverse dynamics using pinocchio::rnea
    trajPtr->compute(z, false);
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