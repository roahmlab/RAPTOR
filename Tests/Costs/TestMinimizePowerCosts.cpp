#define BOOST_TEST_MODULE MinimizePowerTest
#include <boost/test/included/unit_test.hpp>

#include "MinimizePower.h"
#include "Polynomials.h"
#include <chrono>

using namespace RAPTOR;

BOOST_AUTO_TEST_SUITE(MinimizePowerTest)

// test gradient
BOOST_AUTO_TEST_CASE(GradientTest) {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 5, model.nv, TimeDiscretization::Chebyshev, 3);
    std::shared_ptr<InverseDynamics> idPtr_ = std::make_shared<InverseDynamics>(model, trajPtr_);
    std::shared_ptr<Costs> costPtr_ = std::make_shared<MinimizePower>(trajPtr_, idPtr_);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr_->varLength);
    
    costPtr_->compute(z, true, false);

    const Eigen::VectorXd J_analytical = costPtr_->grad_f;
    Eigen::VectorXd J_numerical = Eigen::VectorXd::Zero(J_analytical.size());
    for (int i = 0; i < z.size(); i++) {
        Eigen::VectorXd q_plus = z;
        q_plus(i) += 1e-8;
        costPtr_->compute(q_plus, false);
        const double f_plus = costPtr_->f;
        Eigen::VectorXd q_minus = z;
        q_minus(i) -= 1e-8;
        costPtr_->compute(q_minus, false);
        const double f_minus = costPtr_->f;
        J_numerical(i) = (f_plus - f_minus) / 2e-8;
    }

    BOOST_CHECK_SMALL((J_analytical - J_numerical).norm(), 1e-4);
}

// test hessian
BOOST_AUTO_TEST_CASE(HessianTest) {
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Polynomials>(2.0, 5, model.nv, TimeDiscretization::Uniform, 3);
    std::shared_ptr<InverseDynamics> idPtr_ = std::make_shared<InverseDynamics>(model, trajPtr_);
    std::shared_ptr<Costs> costPtr_ = std::make_shared<MinimizePower>(trajPtr_, idPtr_);

    // compute a valid transform using forward kinematics
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(trajPtr_->varLength);

    costPtr_->compute(z, true, true);

    const Eigen::MatrixXd H_analytical = costPtr_->hess_f;
    for (int i = 0; i < z.size(); i++) {
        Eigen::VectorXd q_plus = z;
        q_plus(i) += 1e-8;
        costPtr_->compute(q_plus, true, false);
        const Eigen::VectorXd J_plus = costPtr_->grad_f;
        Eigen::VectorXd q_minus = z;
        q_minus(i) -= 1e-8;
        costPtr_->compute(q_minus, true, false);
        const Eigen::VectorXd J_minus = costPtr_->grad_f;
        const Eigen::VectorXd H_numerical_row = (J_plus - J_minus) / 2e-8;
        
        BOOST_CHECK_SMALL((H_analytical.row(i) - H_numerical_row.transpose()).norm(), 1e-5);
    }
}

BOOST_AUTO_TEST_SUITE_END()
