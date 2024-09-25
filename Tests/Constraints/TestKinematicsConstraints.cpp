#include "KinematicsConstraints.h"
// #include "Plain.h"
#include "Polynomials.h"
#include <chrono>

using namespace RAPTOR;

int main() {
    // Define robot model
    // const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    pinocchio::ModelTpl<float> model = model_double.cast<float>();

    // std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(model.nv);
    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Chebyshev, 3);
    ForwardKinematicsSolver fkSolver(&model);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXf z = M_2_PI * Eigen::VectorXf::Random(trajPtr_->varLength).array() - M_PI;
    int start = 0;
    int end = model.getJointId("joint_7");
    fkSolver.compute(start, end, z);

    KinematicsConstraints kc(trajPtr_, &model, end, 6, fkSolver.getTransform());

    // simple test when difference is small
    Eigen::VectorXf z_test = z.array() + 1e-6;
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
    const Eigen::MatrixXf J_analytical = kc.pg_pz;
    Eigen::MatrixXf J_numerical = Eigen::MatrixXf::Zero(J_analytical.rows(), J_analytical.cols());
    for (int i = 0; i < z_test.size(); i++) {
        Eigen::VectorXf q_plus = z_test;
        q_plus(i) += 1e-7;
        kc.compute(q_plus, false);
        const Eigen::VectorXf g_plus = kc.g;
        Eigen::VectorXf q_minus = z_test;
        q_minus(i) -= 1e-7;
        kc.compute(q_minus, false);
        const Eigen::VectorXf g_minus = kc.g;
        J_numerical.col(i) = (g_plus - g_minus) / 2e-7;
    }

    // std::cout << "Analytical gradient: " << std::endl << J_analytical << std::endl << std::endl;
    // std::cout << "Numerical gradient: " << std::endl << J_numerical << std::endl << std::endl;
    std::cout << J_analytical - J_numerical << std::endl << std::endl;

    // test hessian
    Eigen::Array<Eigen::MatrixXf, 1, Eigen::Dynamic> H_analytical = kc.pg_pz_pz;
    for (int i = 0; i < z_test.size(); i++) {
        Eigen::VectorXf q_plus = z_test;
        q_plus(i) += 1e-7;
        kc.compute(q_plus, true, false);
        const Eigen::MatrixXf J_plus = kc.pg_pz;
        Eigen::VectorXf q_minus = z_test;
        q_minus(i) -= 1e-7;
        kc.compute(q_minus, true, false);
        const Eigen::MatrixXf J_minus = kc.pg_pz;
        const Eigen::MatrixXf H_numerical_row = (J_plus - J_minus) / 2e-7;

        for (int j = 0; j < 3; j++) {
            std::cout << H_analytical(j).row(i) - H_numerical_row.row(j) << std::endl;
        }
    }

    return 0;
}
