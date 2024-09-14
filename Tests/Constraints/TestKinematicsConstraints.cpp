#define BOOST_TEST_MODULE KinematicsConstraintsTest
#include <boost/test/included/unit_test.hpp>

#include "KinematicsConstraints.h"
// #include "Plain.h"
#include "Polynomials.h"
#include <chrono>

using namespace RAPTOR;


BOOST_AUTO_TEST_SUITE(KinematicsConstraintsTest)
BOOST_AUTO_TEST_CASE(owngradientTest)
{
    // Define robot model
    // const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(model.nv);
    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Chebyshev, 3);
    ForwardKinematicsSolver fkSolver(&model);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = 2 * M_PI * Eigen::VectorXd::Random(trajPtr_->varLength).array() - M_PI;
    int start = 0;
    int end = model.getJointId("joint_7");
    fkSolver.compute(start, end, z);

    KinematicsConstraints kc(trajPtr_, &model, end, 6, fkSolver.getTransform());

    // simple test when difference is small
    Eigen::VectorXd z_test = z.array() + 1e-6;
    // kc.compute(z_test, false);
    // std::cout << kc.g << std::endl << std::endl;

    // simple test when difference is large
    // z_test = z.array() + 1.0;

    // auto start_clock = std::chrono::high_resolution_clock::now();
    // kc.compute(z_test, true, true);
    // auto end_clock = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_clock - start_clock);
    // std::cout << "Time taken (including gradient and hessian): " << duration.count() << " microseconds" << std::endl;

    // std::cout << kc.g << std::endl << std::endl;

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
    BOOST_CHECK_SMALL((J_analytical - J_numerical).norm(), 1e-10);
    }

    
    // test hessian

    BOOST_AUTO_TEST_CASE(ownHessianTest){
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(model.nv);
    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Chebyshev, 3);
    ForwardKinematicsSolver fkSolver(&model);


    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = 2 * M_PI * Eigen::VectorXd::Random(trajPtr_->varLength).array() - M_PI;
    int start = 0;
    int end = model.getJointId("joint_7");
    fkSolver.compute(start, end, z);

    // check error
    bool hasError = false;
    double max_diff = 0.0;
    std::stringstream error_message;

    KinematicsConstraints kc(trajPtr_, &model, end, 6, fkSolver.getTransform());

    // simple test when difference is small
    Eigen::VectorXd z_test = z.array() + 1e-6;
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

        // check error 
        for (int j = 0; j < 3; j++) {
            // sstd::cout << H_analytical(j).row(i) - H_numerical_row.row(j) << std::endl;

            //check each loop
            // BOOST_CHECK_SMALL((H_analytical(j).row(i) - H_numerical_row.row(j) ).norm(), 1e-10);

            //check one time
            double diff = (H_analytical(j).row(i) - H_numerical_row.row(j) ).norm();
            if (diff >1e-10){
                hasError = true;
                if (diff >max_diff) max_diff = diff;
                error_message << "error found at i=" << i << ", j=" << j 
                              << " with difference: " << diff << "\n";

            }
            
        }

        
    }
    // bool hasError = false;
    // double max_diff = 0.0;
    // for (int i = 0; i < z_test.size(); i++) {
    //     for (int j = 0; j < 3; j++){
    //         double diff = (H_analytical(j).row(i) - H_numerical_row.row(j) ).norm()
    //         if (diff >1e-10){
    //             hasError = true;
    //             if (diff >max_diff) max_diff = diff;
    //         }
    //         error_message << "Discrepancy found at i=" << i << ", j=" << j 
    //                           << " with difference: " << diff << "\n";
    //     }
    // }
    BOOST_CHECK_MESSAGE(!hasError, "Hessian discrepancies found:\n" + error_message.str());


  
}
BOOST_AUTO_TEST_SUITE_END()
