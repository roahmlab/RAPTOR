#include "KinovaConstants.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "FixedFrequencyFourierCurves.h"
#include "BezierCurves.h"
#include "Plain.h"

#include "RegressorInverseDynamics.h"

using namespace IDTO;
using namespace Kinova;
// using namespace Ipopt;

int main() {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);
    
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.81;
    model.friction.setZero();
    // model.damping.setZero();
    // model.rotorInertia.setZero();

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    std::shared_ptr<Trajectories> traj = std::make_shared<FixedFrequencyFourierCurves>(1.0, 
                                                                                       5, 
                                                                                       model.nv, 
                                                                                       TimeDiscretization::Uniform, 
                                                                                       4);
    // std::shared_ptr<Trajectories> traj = std::make_shared<BezierCurves>(1.0, 
    //                                                                     5, 
    //                                                                     model.nv, 
    //                                                                     TimeDiscretization::Uniform, 
    //                                                                     2);
    // std::shared_ptr<Trajectories> traj = std::make_shared<Plain>(model.nv);

    RegressorInverseDynamics rid(model, jtype, traj);
    InverseDynamics id(model, traj);

    Eigen::VectorXd z = Eigen::VectorXd::Random(traj->varLength);

    auto start = std::chrono::high_resolution_clock::now();
    rid.compute(z, true);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Regressor Inverse Dynamics: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    id.compute(z, false);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Inverse Dynamics: " << duration.count() << " microseconds" << std::endl;

    for (int i = 0; i < rid.tau.size(); i++) {
        std::cout << rid.tau(i).transpose() - id.tau(i).transpose() << std::endl;
    }

    // compute numerical gradient by central finite difference
    for (int j = 0; j < traj->N; j++) {
        for (int i = 0; i < z.size(); i++) {
            Eigen::VectorXd z_plus = z;
            Eigen::VectorXd z_minus = z;
            z_plus(i) += 1e-8;
            z_minus(i) -= 1e-8;

            rid.compute(z_plus, false);
            Eigen::VectorXd tau_plus = rid.tau(j);

            rid.compute(z_minus, false);
            Eigen::VectorXd tau_minus = rid.tau(j);

            Eigen::VectorXd grad = (tau_plus - tau_minus) / 2e-8;
            std::cout << (grad - rid.ptau_pz(j).col(i)).norm() << std::endl;
            // std::cout <<  << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}