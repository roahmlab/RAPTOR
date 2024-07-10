#include "KinematicsConstraints.h"
// #include "Plain.h"
#include "Polynomials.h"

#include <chrono>

using namespace IDTO;

int main() {
    // Define robot model
    // const std::string urdf_filename = "../Examples/Digit/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    // jtype << 4, 5, 6, 1, 2, 3, 
    //          3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3,
    //          3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3;
    jtype << 3, 3, 3, 3, 3, 3, 3;

    // std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(model.nv);
    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Chebyshev, 3);
    ForwardKinematicsSolver fkSolver(&model, jtype);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = 2 * M_PI * Eigen::VectorXd::Random(trajPtr_->varLength).array() - M_PI;
    int start = 0;
    int end = model.getJointId("joint_7");
    fkSolver.compute(start, end, z);

    KinematicsConstraints kc(trajPtr_, &model, jtype, end, 6, fkSolver.getTransform());

    // simple test when difference is small
    Eigen::VectorXd z_test = z.array() + 1e-6;
    kc.compute(z_test, false);
    std::cout << kc.g << std::endl << std::endl;

    // simple test when difference is large
    z_test = z.array() + 1.0;
    auto start_clock = std::chrono::high_resolution_clock::now();
    kc.compute(z_test, true, true);
    auto end_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_clock - start_clock);
    std::cout << "Time taken (including gradient and hessian): " << duration.count() << " microseconds" << std::endl;

    std::cout << kc.g << std::endl << std::endl;

    // test gradient
    const Eigen::MatrixXd J_analytical = kc.pg_pz;
    Eigen::MatrixXd J_numerical = Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());
    for (int i = 0; i < z_test.size(); i++) {
        Eigen::VectorXd q_plus = z_test;
        q_plus(i) += 1e-7;
        kc.compute(q_plus, false);
        const Eigen::VectorXd g_plus = kc.g;
        Eigen::VectorXd q_minus = z_test;
        q_minus(i) -= 1e-7;
        kc.compute(q_minus, false);
        const Eigen::VectorXd g_minus = kc.g;
        J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    }

    // std::cout << "Analytical gradient: " << std::endl << J_analytical << std::endl << std::endl;
    // std::cout << "Numerical gradient: " << std::endl << J_numerical << std::endl << std::endl;
    std::cout << J_analytical - J_numerical << std::endl << std::endl;

    // test hessian
    Eigen::Array<Eigen::MatrixXd, 1, Eigen::Dynamic> H_analytical = kc.pg_pz_pz;
    for (int i = 0; i < z_test.size(); i++) {
        Eigen::VectorXd q_plus = z_test;
        q_plus(i) += 1e-7;
        kc.compute(q_plus, true, false);
        const Eigen::MatrixXd J_plus = kc.pg_pz;
        Eigen::VectorXd q_minus = z_test;
        q_minus(i) -= 1e-7;
        kc.compute(q_minus, true, false);
        const Eigen::MatrixXd J_minus = kc.pg_pz;
        const Eigen::MatrixXd H_numerical_row = (J_plus - J_minus) / 2e-7;

        for (int j = 0; j < 3; j++) {
            std::cout << H_analytical(j).row(i) - H_numerical_row.row(j) << std::endl;
        }
    }

    return 0;
}