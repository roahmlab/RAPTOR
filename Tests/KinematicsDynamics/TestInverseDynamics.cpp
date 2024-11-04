#define BOOST_TEST_MODULE InverseDynamicsTest
#include <boost/test/included/unit_test.hpp>

#include "InverseDynamics.h"
#include "Polynomials.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(InverseDynamicsTest)

// test gradient
BOOST_AUTO_TEST_CASE(GradientTest) {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Chebyshev, 3);
    std::shared_ptr<InverseDynamics> idPtr_ = std::make_shared<InverseDynamics>(model, trajPtr_);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr_->varLength);
    
    idPtr_->compute(z, true, false);

    for (int i = 0; i < idPtr_->N; i++) {
        const Eigen::MatrixXd J_analytical = idPtr_->ptau_pz(i);
        Eigen::MatrixXd J_numerical = Eigen::MatrixXd::Zero(J_analytical.rows(), J_analytical.cols());
        for (int j = 0; j < z.size(); j++) {
            Eigen::VectorXd q_plus = z;
            q_plus(j) += 1e-8;
            idPtr_->compute(q_plus, false);
            const Eigen::VectorXd f_plus = idPtr_->tau(i);
            Eigen::VectorXd q_minus = z;
            q_minus(j) -= 1e-8;
            idPtr_->compute(q_minus, false);
            const Eigen::VectorXd f_minus = idPtr_->tau(i);
            J_numerical.col(j) = (f_plus - f_minus) / 2e-8;
        }

        BOOST_CHECK_SMALL((J_analytical - J_numerical).norm(), 1e-5);
    }
}

// test hessian
BOOST_AUTO_TEST_CASE(HessianTest) {
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 10, model.nv, TimeDiscretization::Uniform, 3);
    std::shared_ptr<InverseDynamics> idPtr_ = std::make_shared<InverseDynamics>(model, trajPtr_);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr_->varLength);

    idPtr_->compute(z, true, true);

    for (int i = 0; i < idPtr_->N; i++) {
        const Eigen::MatrixXd H_analytical = idPtr_->ptau_pz_pz(6, i);
        for (int j = 0; j < z.size(); j++) {
            Eigen::VectorXd q_plus = z;
            q_plus(j) += 1e-8;
            idPtr_->compute(q_plus, true, false);
            const Eigen::MatrixXd J_plus = idPtr_->ptau_pz(i);
            Eigen::VectorXd q_minus = z;
            q_minus(j) -= 1e-8;
            idPtr_->compute(q_minus, true, false);
            const Eigen::MatrixXd J_minus = idPtr_->ptau_pz(i);
            const Eigen::VectorXd H_numerical_row = (J_plus.row(6) - J_minus.row(6)) / 2e-8;
            BOOST_CHECK_SMALL( (H_analytical.row(j) - H_numerical_row.transpose()).norm(), 1e-5);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
