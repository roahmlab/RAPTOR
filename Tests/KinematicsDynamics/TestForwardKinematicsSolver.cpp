#include "ForwardKinematics.h"

using namespace IDTO;

int main() {
    // Define robot model
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3,
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3;

    ForwardKinematicsSolver fkSolver(&model, jtype);

    // set joint angles
    std::srand(std::time(nullptr));
    Eigen::VectorXd q = 2 * M_PI * Eigen::VectorXd::Random(model.nq).array() - M_PI;

    // compute forward kinematics using pinocchio
    auto start_clock = std::chrono::high_resolution_clock::now();
    pinocchio::forwardKinematics(model, data, q);
    auto stop_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "Pinocchio FK: " << duration.count() << " nanoseconds" << std::endl;

    // set the start and end joint
    int start = 0;
    int end = model.getJointId("left_toe_B");

    // compute forward kinematics using IDTO
    start_clock = std::chrono::high_resolution_clock::now();
    fkSolver.compute(start, end, q);
    stop_clock = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "IDTO FK: " << duration.count() << " nanoseconds" << std::endl;

    // compare the results
    std::cout << "Pinocchio: " << data.oMi[model.getJointId("left_toe_B")].translation().transpose() << std::endl;
    std::cout << "IDTO: " << fkSolver.getTranslation().transpose() << std::endl;

    return 0;
}