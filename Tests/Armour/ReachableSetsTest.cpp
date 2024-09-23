#include "ReachableSets.h"
#include "ArmourBezierCurves.h"
#include "InverseDynamics.h"

using namespace RAPTOR;
using namespace Armour;

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

// INITIALIZATION
    // read robot model and info
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithoutGripperInfo.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr_ = 
        std::make_shared<RobotInfo>(robot_model_file, robot_info_file);

    // turn off tracking error for validation
    robotInfoPtr_->ultimate_bound_info.eps = 0;
    robotInfoPtr_->ultimate_bound_info.qe = 0;
    robotInfoPtr_->ultimate_bound_info.qde = 0;
    robotInfoPtr_->ultimate_bound_info.qdae = 0;
    robotInfoPtr_->ultimate_bound_info.qddae = 0;

    // turn off friction for validation
    robotInfoPtr_->model.friction.setZero();
    
    // create a trajectory instance (compute trajectory on continuous time intervals)
        // initial conditions of the trajectory
    const Eigen::VectorXf q0 = Eigen::VectorXf::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXf q_d0 = Eigen::VectorXf::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXf q_dd0 = Eigen::VectorXf::Random(robotInfoPtr_->num_motors);

        // trajectory parameters and their ranges
    const Eigen::VectorXf k_center = Eigen::VectorXf::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXf k_range = M_PI / 24 * Eigen::VectorXf::Ones(robotInfoPtr_->num_motors);

        // trajectory duration
    const float duration = 2.0;

    std::shared_ptr<BezierCurveInterval> trajPtr_ = 
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_->ultimate_bound_info);
    
    // create a KinematicsDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<KinematicsDynamics> kdPtr = 
        std::make_shared<KinematicsDynamics>(robotInfoPtr_, trajPtr_);

// COMPUTATION IN ARMOUR 
    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    GenerateJRS(robotInfoPtr_, trajPtr_);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate JRS: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;
    
    // generate Link and Torque PZs
    auto start2 = std::chrono::high_resolution_clock::now();
    GenerateLinkAndTorquePZs(robotInfoPtr_, trajPtr_, kdPtr);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate Link and Torque PZs: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() 
              << " ms" << std::endl;

    // // compute Robust Input Bounds
    // auto start3 = std::chrono::high_resolution_clock::now();
    // ComputeRobustInputBounds(robotInfoPtr_, trajPtr_, kdPtr);
    // auto end3 = std::chrono::high_resolution_clock::now();
    // std::cout << "Time taken to compute Robust Input Bounds: " 
    //           << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count() 
    //           << " ms" << std::endl;

// VALIDATION
    // randomly choose a trajectory parameter inside the range
    const Eigen::VectorXf factor = Eigen::VectorXf::Random(robotInfoPtr_->num_motors); // [-1, 1]
    const Eigen::VectorXf k = k_center + k_range.cwiseProduct(factor);

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
            (i + static_cast<float>(std::rand()) / RAND_MAX) * 
                (duration / trajPtr_->num_time_steps);
    }

    trajDiscretePtr_->compute(k, false);

    // validate the joint trajectory reachable sets
        // cos of q (joint trajectory)
    for (int i = 0; i < trajPtr_->num_time_steps; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            // slice the JRS PZ at k (the actual trajectory parameter)
            const Interval cos_q_range = trajPtr_->cos_q_des(j, i).slice(factor);
            const float actual_cos_q = cosf(trajDiscretePtr_->q(i)(j));

            // check if actual_cos_q is inside the range
            if (actual_cos_q < cos_q_range.lower() - 1e-5 || 
                actual_cos_q > cos_q_range.upper() + 1e-5) {
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
            const Interval sin_q_range = trajPtr_->sin_q_des(j, i).slice(factor);
            const float actual_sin_q = sinf(trajDiscretePtr_->q(i)(j));
            if (actual_sin_q < sin_q_range.lower() - 1e-5 || 
                actual_sin_q > sin_q_range.upper() + 1e-5) {
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
            const float actual_q_d = trajDiscretePtr_->q_d(i)(j);
            if (actual_q_d < q_d_range.lower() - 1e-4 || 
                actual_q_d > q_d_range.upper() + 1e-4) {
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
            const float actual_q_dd = trajDiscretePtr_->q_dd(i)(j);
            if (actual_q_dd < q_dd_range.lower() - 1e-3 || 
                actual_q_dd > q_dd_range.upper() + 1e-3) {
                std::cerr << "Validation failed for q_dd at time step " << i 
                          << " for motor " << j << ": "
                          << actual_q_dd << " not in [ " 
                          << q_dd_range.lower() << ", " 
                          << q_dd_range.upper() << " ]" << std::endl;
            }
        }
    }

    // validate the torque PZs
    std::shared_ptr<InverseDynamics> idPtr_ = 
        std::make_shared<InverseDynamics>(robotInfoPtr_->model, trajDiscretePtr_);
    idPtr_->compute(k, false);
    
    for (int i = 0; i < idPtr_->N; i++) {
        for (int j = 0; j < robotInfoPtr_->num_motors; j++) {
            const Interval torqueRange = kdPtr->torque_nom(j, i).slice(factor);
            const float actualTorque = idPtr_->tau(i)(j);

            if (actualTorque < torqueRange.lower() - 1e-3 || 
                actualTorque > torqueRange.upper() + 1e-3) {
                std::cerr << "Validation failed for tau at time step " << i 
                          << " for motor " << j << ": "
                          << actualTorque << " not in [ " 
                          << torqueRange.lower() << ", " 
                          << torqueRange.upper() << " ]" << std::endl;
            }
        }
    }

    return 0;
}