#include "KinematicsConstraints.h"
#include "Plain.h"

using namespace IDTO;

int main() {
    // Define robot model
    const std::string urdf_filename = "../Examples/Digit/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
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

    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(model.nv);
    ForwardKinematicsSolver fkSolver(&model, jtype);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd q = 2 * M_PI * Eigen::VectorXd::Random(model.nq).array() - M_PI;
    int start = 0;
    int end = model.getJointId("left_toe_B");
    fkSolver.compute(start, end, q);

    KinematicsConstraints kc(trajPtr_, &model, jtype, end, 0, fkSolver.getTransform());

    Eigen::VectorXd q_test = q.array() + 0.1;
    kc.compute(q_test, true);

    std::cout << kc.g << std::endl << std::endl;

    // test gradient
    Eigen::MatrixXd J_analytical = kc.pg_pz;
    Eigen::MatrixXd J_numerical = Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());
    for (int i = 0; i < q_test.size(); i++) {
        Eigen::VectorXd q_plus = q_test;
        q_plus(i) += 1e-7;
        kc.compute(q_plus, false);
        const Eigen::VectorXd g_plus = kc.g;
        Eigen::VectorXd q_minus = q_test;
        q_minus(i) -= 1e-7;
        kc.compute(q_minus, false);
        const Eigen::VectorXd g_minus = kc.g;
        J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    }

    // std::cout << "Analytical gradient: " << std::endl << J_analytical << std::endl << std::endl;
    // std::cout << "Numerical gradient: " << std::endl << J_numerical << std::endl << std::endl;
    std::cout << J_analytical - J_numerical << std::endl;

    return 0;
}