#include "DigitMultipleStepOptimizer.h"

using namespace RAPTOR;
using namespace Digit;
using namespace Ipopt;

const std::string filepath = "../Examples/Digit/data/";

int main() {
    const std::string urdf_filename = "../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;
    
    // ignore friction for now
    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero();
    
    const int numSteps = 2;
    const int degree = 5;
    const double T = 0.4;
    const double FPS = 60.0;
    const int N = T * FPS;
    GaitParameters gp;
    
    const Eigen::VectorXd solution = Utils::initializeEigenMatrixFromFile(filepath + "solution_upstairs.txt");

    std::ofstream trajectories(filepath + "full-trajectories-upstairs.txt");
    trajectories << std::setprecision(20);

    // setup optimizers
    std::vector<SmartPtr<DigitSingleStepOptimizer>> testnlps;
    testnlps.reserve(numSteps);

    int offset = 0;
    Transform previousStandingFootTransform;
    
    for (int step = 0; step < numSteps; step++) {
        testnlps.push_back(new DigitSingleStepOptimizer());

        // Eigen::VectorXd z((degree + 1) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS);
        // for (int i = 0; i < z.size(); i++) {
        //     z(i) = solution(offset + i);
        // }
        // offset += z.size();

        char stanceLeg = (step % 2 == 0) ? 'L' : 'R';

        Eigen::VectorXd z;
        if (stanceLeg == 'R') {
            z = switchSolutionFromLeftToRight(solution, degree);
        }
        else {
            z = solution;
        }

        Transform stanceFootTransform(3, -M_PI / 2);
        if (step > 0) {
            stanceFootTransform = previousStandingFootTransform;
        }

        auto& testnlp = testnlps[step];
        testnlp->set_parameters(z,
                                T,
                                N,
                                TimeDiscretization::Uniform,
                                degree,
                                model,
                                gp,
                                stanceLeg,
                                stanceFootTransform);  

        // evaluate everything
        Index n, m, nnz_jac_g, nnz_h_lag;
        TNLP::IndexStyleEnum index_style;
        testnlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
        Number ztry[testnlp->numVars], x_l[testnlp->numVars], x_u[testnlp->numVars];
        Number g[testnlp->numCons], g_lb[testnlp->numCons], g_ub[testnlp->numCons];
        for (int i = 0; i < testnlp->numVars; i++) {
            ztry[i] = z(i);
        }
        testnlp->get_bounds_info(testnlp->numVars, x_l, x_u, testnlp->numCons, g_lb, g_ub);
        testnlp->eval_g(testnlp->numVars, ztry, false, testnlp->numCons, g);
        testnlp->summarize_constraints(testnlp->numCons, g, false);

        // print trajectories to file
        for (int j = 0; j < testnlp->cidPtr_->N; j++) {
            trajectories << testnlp->cidPtr_->q(j).transpose() << std::endl;
        }

        // compute swing foot end config, which is stance foot start config for next step
        Eigen::VectorXd swingfoot_xyzrpy;
        for (size_t j = 0; j < testnlp->constraintsNameVec_.size(); j++) {
            if (testnlp->constraintsNameVec_[j] == "customized constraints") {
                const auto& constraintsPtr_ = testnlp->constraintsPtrVec_[j];
                const auto& customizedConstraintsPtr = dynamic_cast<DigitCustomizedConstraints*>(constraintsPtr_.get());

                swingfoot_xyzrpy = 
                    customizedConstraintsPtr->swingfoot_xyzrpy.col(customizedConstraintsPtr->swingfoot_xyzrpy.cols() - 1);
                
                break;
            }
        }

        previousStandingFootTransform = Transform(
            Eigen::Vector3d(swingfoot_xyzrpy.tail(3)), 
            Eigen::Vector3d(swingfoot_xyzrpy.head(3)));
    }
}