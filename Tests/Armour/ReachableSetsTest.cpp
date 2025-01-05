#include "PZDynamics.h"
#include "ArmourBezierCurves.h"
#include "InverseDynamics.h"
#include "CustomizedInverseDynamics.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

// INITIALIZATION
    // read robot model and info
        // Note that you need to comment line 166 for the following robot test
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithoutGripperInfo.yaml";

        // Note that you need to uncomment line 165 for the following robot test
    // const std::string robot_model_file = "../Robots/kinova-gen3/kinova_grasp.urdf";
    // const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaSuctionCup.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr_ = 
        std::make_shared<RobotInfo>(robot_model_file, robot_info_file);

    // turn off friction for validation
    robotInfoPtr_->model.friction.setZero();

    // turn off tracking error for validation
    robotInfoPtr_->ultimate_bound_info.eps = 0;
    robotInfoPtr_->ultimate_bound_info.qe = 0;
    robotInfoPtr_->ultimate_bound_info.qde = 0;
    robotInfoPtr_->ultimate_bound_info.qdae = 0;
    robotInfoPtr_->ultimate_bound_info.qddae = 0;
    
    // create a trajectory instance (compute trajectory on continuous time intervals)
        // initial conditions of the trajectory
    const Eigen::VectorXd q0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_d0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_dd0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);

        // trajectory parameters and their ranges
    const Eigen::VectorXd k_center = Eigen::VectorXd::Zero(robotInfoPtr_->num_motors);
    const Eigen::VectorXd k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr_->num_motors);

        // trajectory duration
    const double duration = 2.0;
    const size_t num_time_steps = 256;

    std::shared_ptr<BezierCurveInterval> trajPtr_ = 
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_,
            num_time_steps);

    // create a PZDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<PZDynamics> dynPtr_ = 
        std::make_shared<PZDynamics>(robotInfoPtr_, trajPtr_);

// COMPUTATION IN ARMOUR 
    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;

// VALIDATION
    // randomly choose a trajectory parameter inside the range
    const Eigen::VectorXd factor = Eigen::VectorXd::Zero(robotInfoPtr_->num_motors); // [-1, 1]
    const Eigen::VectorXd k = k_center + k_range.cwiseProduct(factor);

    // create a trajectory instance (compute trajectory on discrete time instances)
    ArmourTrajectoryParameters atp;
    atp.q0 = q0;
    atp.q_d0 = q_d0;
    atp.q_dd0 = q_dd0;

    std::shared_ptr<Trajectories> trajDiscretePtr_ = 
        std::make_unique<ArmourBezierCurves>(
            duration, trajPtr_->num_time_steps, robotInfoPtr_->num_motors, Uniform, atp);

    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        // randomly sample inside each time interval
        trajDiscretePtr_->tspan(i) = 
            (i + static_cast<double>(std::rand()) / RAND_MAX) * 
                (duration / trajPtr_->num_time_steps);
    }

    // validate the joint trajectory reachable sets
    trajDiscretePtr_->compute(k, false);

        // cos of q (joint trajectory)
    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            const PZSparse test_pz = cos(trajPtr_->q_des(j, i));
            const Interval cos_q_range = test_pz.slice(factor);
            const double actual_cos_q = std::cos(trajDiscretePtr_->q(i)(j));

            // check if actual_cos_q is inside the range
            if (actual_cos_q < cos_q_range.lower() || 
                actual_cos_q > cos_q_range.upper()) {
                std::cerr << "Validation failed for cos(q) at time step " << i 
                          << " for motor " << j << ": "
                          << actual_cos_q << " not in [ " 
                          << cos_q_range.lower() << ", " 
                          << cos_q_range.upper() << " ]" << std::endl;
            }
        }
    }

        // sin of q (joint trajectory)
    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            const PZSparse test_pz = sin(trajPtr_->q_des(j, i));
            const Interval sin_q_range = test_pz.slice(factor);
            const double actual_sin_q = std::sin(trajDiscretePtr_->q(i)(j));

            if (actual_sin_q < sin_q_range.lower() || 
                actual_sin_q > sin_q_range.upper()) {
                std::cerr << "Validation failed for sin(q) at time step " << i 
                          << " for motor " << j << ": "
                          << actual_sin_q << " not in [ " 
                          << sin_q_range.lower() << ", " 
                          << sin_q_range.upper() << " ]" << std::endl;
            }
        }
    }

        // q_d (joint velocity)
    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            const Interval q_d_range = trajPtr_->qd_des(j, i).slice(factor);
            const double actual_q_d = trajDiscretePtr_->q_d(i)(j);
            if (actual_q_d < q_d_range.lower() || 
                actual_q_d > q_d_range.upper()) {
                std::cerr << "Validation failed for q_d at time step " << i 
                          << " for motor " << j << ": "
                          << actual_q_d << " not in [ " 
                          << q_d_range.lower() << ", " 
                          << q_d_range.upper() << " ]" << std::endl;
            }
        }
    }

        // q_dd (joint acceleration)
    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            Interval q_dd_range = trajPtr_->qdda_des(j, i).slice(factor);
            const double actual_q_dd = trajDiscretePtr_->q_dd(i)(j);
            if (actual_q_dd < q_dd_range.lower() || 
                actual_q_dd > q_dd_range.upper()) {
                std::cerr << "Validation failed for q_dd at time step " << i 
                          << " for motor " << j << ": "
                          << actual_q_dd << " not in [ " 
                          << q_dd_range.lower() << ", " 
                          << q_dd_range.upper() << " ]" << std::endl;
            }
        }
    }

    // validate the torque PZs
    Eigen::VectorXi jtype = convertPinocchioJointType(robotInfoPtr_->model);

    // Note: uncomment the following line if using KinovaSuctionCup.yaml
    // jtype(jtype.size() - 1) = 0;

    std::shared_ptr<CustomizedInverseDynamics> cidPtr_ = 
        std::make_shared<CustomizedInverseDynamics>(robotInfoPtr_->model, trajDiscretePtr_, jtype);
    cidPtr_->compute(k, false);
    
    for (int i = 0; i < cidPtr_->N; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            const Interval torqueRange = dynPtr_->data_sparses[i].tau(j).slice(factor);
            const double actualTorque = cidPtr_->tau(i)(j);

            if (actualTorque < torqueRange.lower() || 
                actualTorque > torqueRange.upper()) {
                std::cerr << "Validation failed for tau at time step " << i 
                          << " for motor " << j << ": "
                          << actualTorque << " not in [ " 
                          << torqueRange.lower() << ", " 
                          << torqueRange.upper() << " ]" << std::endl;
            }
        }
    }

    // validate the contact force PZs
    if (robotInfoPtr_->num_joints > robotInfoPtr_->num_motors) {
        for (int i = 0; i < cidPtr_->N; i++) {
            const Eigen::Vector<double, 6>& actualLambda = cidPtr_->lambda(i);
            for (int j = 0; j < 3; j++) {
                const Interval forceRange = dynPtr_->data_sparses[i].f[dynPtr_->model_sparses[i].nv].linear()(j).slice(factor);
                const Interval momentRange = dynPtr_->data_sparses[i].f[dynPtr_->model_sparses[i].nv].angular()(j).slice(factor);            
                if (actualLambda(j + 3) < forceRange.lower() || 
                    actualLambda(j + 3) > forceRange.upper()) {
                    std::cerr << "Validation failed for contact force at time step " << i 
                            << " for direction " << j << ": "
                            << actualLambda(j+3) << " not in [ " 
                            << forceRange.lower() << ", " 
                            << forceRange.upper() << " ]" << std::endl;
                }
                if (actualLambda(j) < momentRange.lower() || 
                    actualLambda(j) > momentRange.upper()) {
                    std::cerr << "Validation failed for contact moment at time step " << i 
                            << " for direction " << j << ": "
                            << actualLambda(j) << " not in [ " 
                            << momentRange.lower() << ", " 
                            << momentRange.upper() << " ]" << std::endl;
                }
            }
        }
    }

    // // validate that torque_int overapproximate torque_nom
    // // you need to comment out line 44 - 47 in ReachableSets.cpp
    // for (int i = 0; i < cidPtr_->N; i++) {
    //     for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
    //         const Interval torqueRange = dynPtr_->torque_int(j, i).slice(factor);
    //         const Interval torqueNomRange = dynPtr_->torque_nom(j, i).slice(factor);
    //         if (!((torqueRange.lower() < torqueNomRange.lower() - 1e-3) &&
    //               (torqueRange.upper() > torqueNomRange.upper() + 1e-3)))
    //         {
    //             std::cerr << "Validation failed for torque_int at time step " << i
    //                       << " for motor " << j << ": "
    //                       << torqueNomRange << " not in"
    //                       << torqueRange << std::endl;
    //         }
    //     }
    // }

    return 0;
}