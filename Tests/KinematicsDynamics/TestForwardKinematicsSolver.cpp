#include "ForwardKinematics.h"
#include <chrono>

using namespace RAPTOR;

int main() {
    // Define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    pinocchio::ModelTpl<float> model = model_double.cast<float>();
    pinocchio::DataTpl<float> data(model);

    ForwardKinematicsSolver fkSolver(&model);

    // set joint angles
    std::srand(std::time(nullptr));
    Eigen::VectorXf q = 2 * M_PI * Eigen::VectorXf::Random(model.nq).array() - M_PI;

    // compute forward kinematics using pinocchio
    auto start_clock = std::chrono::high_resolution_clock::now();
    pinocchio::forwardKinematics(model, data, q);
    auto stop_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "Pinocchio FK: " << duration.count() << " nanoseconds" << std::endl;

    // set the start and end joint
    int start = 0;
    int end = model.getJointId("left_toe_B");

    // compute forward kinematics using RAPTOR
    start_clock = std::chrono::high_resolution_clock::now();
    fkSolver.compute(start, end, q);
    stop_clock = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "RAPTOR FK: " << duration.count() << " nanoseconds" << std::endl;

    // compare the results
    std::cout << "Pinocchio: " << data.oMi[model.getJointId("left_toe_B")].translation().transpose() << std::endl;
    std::cout << "RAPTOR: " << fkSolver.getTranslation().transpose() << std::endl;

    return 0;
}