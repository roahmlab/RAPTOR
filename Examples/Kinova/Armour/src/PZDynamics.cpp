#include "PZDynamics.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

PZDynamics::PZDynamics(const std::shared_ptr<RobotInfo>& robotInfoPtr_input,
			           const std::shared_ptr<BezierCurveInterval>& trajPtr_input) :
    robotInfoPtr_(robotInfoPtr_input),
    trajPtr_(trajPtr_input) {
    // Initialize the sparse model and data
    model_sparses.resize(trajPtr_->num_time_steps);
    data_sparses.resize(trajPtr_->num_time_steps);

    model_sparses_interval.resize(trajPtr_->num_time_steps);
    data_sparses_interval.resize(trajPtr_->num_time_steps);

    for (size_t i = 0; i < trajPtr_->num_time_steps; ++i) {
        model_sparses[i] = robotInfoPtr_->model.cast<PZSparse>();
        data_sparses[i] = pinocchio::DataTpl<PZSparse>(model_sparses[i]);

        model_sparses_interval[i] = robotInfoPtr_->model.cast<PZSparse>();
        data_sparses_interval[i] = pinocchio::DataTpl<PZSparse>(model_sparses_interval[i]);
    }

    torque_radii = Eigen::MatrixXd::Zero(robotInfoPtr_->num_motors, trajPtr_->num_time_steps);
    sphere_radii = Eigen::MatrixXd::Zero(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);

    sample_eigenvalues();
}

void PZDynamics::reset_trajecotry(const std::shared_ptr<BezierCurveInterval>& trajPtr_input) {
    if (trajPtr_input->num_time_steps != trajPtr_->num_time_steps) {
        throw std::invalid_argument("The number of time steps in the trajectory is different from the previous one.");
    }

    trajPtr_ = trajPtr_input;
    sample_eigenvalues();
}

void PZDynamics::sample_eigenvalues(size_t num_samples) {
    pinocchio::Model model = robotInfoPtr_->model;

    if (robotInfoPtr_->num_joints > robotInfoPtr_->num_motors) {
        pinocchio::Model model_reduced;
        std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id = {(pinocchio::JointIndex)model.nv};
        pinocchio::buildReducedModel(model, list_of_joints_to_lock_by_id, Eigen::VectorXd::Zero(model.nv), model_reduced);
        model_reduced.armature = model.armature.head(robotInfoPtr_->num_motors);
        model = model_reduced;
    }
    pinocchio::Data data(model);

    VecX q = VecX::Zero(robotInfoPtr_->num_motors);
    VecX qd = VecX::Zero(robotInfoPtr_->num_motors);
    VecX qdd = VecX::Zero(robotInfoPtr_->num_motors);
    
    double M_max = 0.0;
    double M_min = std::numeric_limits<double>::max();

    std::cout << "Sampling " << num_samples << " configurations..." << std::endl;
    for (int h = 0; h < num_samples; h++) {
        // generate random configurations along the trajectory
        double t = std::rand() / (double)RAND_MAX * trajPtr_->duration;
        VecX k = VecX::Random(robotInfoPtr_->num_motors);
        trajPtr_->computeTrajectories(k, t, q, qd, qdd);
        q += Eigen::VectorXd::Random(robotInfoPtr_->num_motors) * 0.05;

        // compute the inertia mass matrix
        MatX M = pinocchio::crba(model, data, q);
        for (int i = 0; i < M.rows(); i++) {
            for (int j = i + 1; j < M.cols(); j++) {
                M(j, i) = M(i, j);
            }
        }

        // compute the eigenvalues
        Eigen::EigenSolver<MatX> es(M);
        Eigen::VectorXd eigenvalues = es.eigenvalues().real();

        // update the maximum and minimum eigenvalues
        for (int i = 0; i < eigenvalues.size(); i++) {
            M_max = std::max(M_max, eigenvalues(i));
            M_min = std::min(M_min, eigenvalues(i));
        }
    }

    std::cout << "M_max: " << M_max << std::endl;
    std::cout << "M_min: " << M_min << std::endl;

    auto& ultimate_bound_info = robotInfoPtr_->ultimate_bound_info;
    ultimate_bound_info.M_max = M_max;
    ultimate_bound_info.M_min = M_min;

    // recompute the ultimate bound
    robotInfoPtr_->ultimate_bound_info.eps = std::sqrt(2 * ultimate_bound_info.V_m / ultimate_bound_info.M_min);
    ultimate_bound_info.qe = ultimate_bound_info.eps / ultimate_bound_info.Kr;
    ultimate_bound_info.qde = 2 * ultimate_bound_info.eps;
    ultimate_bound_info.qdae = ultimate_bound_info.eps;
    ultimate_bound_info.qddae = 2 * ultimate_bound_info.Kr * ultimate_bound_info.eps;

    std::cout << "Ultimate bound information updated:" << std::endl;
    std::cout << "eps: " << ultimate_bound_info.eps << std::endl;
    std::cout << "qe: " << ultimate_bound_info.qe << std::endl;
    std::cout << "qde: " << ultimate_bound_info.qde << std::endl;
    std::cout << "qdae: " << ultimate_bound_info.qdae << std::endl;
    std::cout << "qddae: " << ultimate_bound_info.qddae << std::endl;
}

void PZDynamics::compute() {
    int t_ind = 0; // openmp loop index
    const size_t num_time_steps = trajPtr_->num_time_steps;
    const size_t num_joints = robotInfoPtr_->num_joints;
    const size_t num_motors = robotInfoPtr_->num_motors;

    // generate joint trajectory reachable sets
    try {
        #pragma omp parallel for shared(trajPtr_) private(t_ind) schedule(static, trajPtr_->num_time_steps / NUM_THREADS)
        for(t_ind = 0; t_ind < num_time_steps; t_ind++) {
            trajPtr_->makePolyZono(t_ind);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error creating JRS! " << e.what() << std::endl;
        throw std::runtime_error("Error creating JRS! Check previous error message!");
    }       

    // generate link and torque PZs
    try {
       #pragma omp parallel for shared(model_sparses, data_sparses, model_sparses_interval, data_sparses_interval, trajPtr_) private(t_ind) schedule(dynamic)
        for(t_ind = 0; t_ind < num_time_steps; t_ind++) {
            Eigen::Vector<PZSparse, Eigen::Dynamic> q(num_joints);
            Eigen::Vector<PZSparse, Eigen::Dynamic> qd(num_joints);
            Eigen::Vector<PZSparse, Eigen::Dynamic> qdd(num_joints);
            for (int i = 0; i < num_joints; i++) {
                if (i < num_motors) {
                    q(i) = trajPtr_->q_des(i, t_ind);   
                    qd(i) = trajPtr_->qda_des(i, t_ind);
                    qdd(i) = trajPtr_->qdda_des(i, t_ind);
                }
                else {
                    q(i) = PZSparse(0);
                    qd(i) = PZSparse(0);
                    qdd(i) = PZSparse(0);
                }
            }

            pinocchio::forwardKinematics(model_sparses[t_ind], data_sparses[t_ind], q);
            pinocchio::updateFramePlacements(model_sparses[t_ind], data_sparses[t_ind]);
            pinocchio::rnea(model_sparses[t_ind], data_sparses[t_ind], q, qd, qdd);
            pinocchio::rnea(model_sparses_interval[t_ind], data_sparses_interval[t_ind], q, qd, qdd);

            for (int i = 0; i < num_motors; i++) {
                data_sparses[t_ind].tau(i) += model_sparses[t_ind].damping(i) * qd(i);
                data_sparses[t_ind].tau(i).reduce();
                data_sparses_interval[t_ind].tau(i) += model_sparses_interval[t_ind].damping(i) * qd(i);
                data_sparses_interval[t_ind].tau(i).reduce();
            }

            for (int i = 0; i < 3; i++) {
                data_sparses[t_ind].f[model_sparses[t_ind].nv].linear()(i).reduce();
                data_sparses[t_ind].f[model_sparses[t_ind].nv].angular()(i).reduce();
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error computing link PZs and nominal torque PZs! " << e.what() << std::endl;
        throw std::runtime_error("Error computing link PZs and nominal torque PZs! Check previous error message!");
    }

    // compute robust input bounds
    const double robust_input_bound = 
        0.5 * robotInfoPtr_->ultimate_bound_info.alpha * 
            (robotInfoPtr_->ultimate_bound_info.M_max - robotInfoPtr_->ultimate_bound_info.M_min) * 
                robotInfoPtr_->ultimate_bound_info.eps;

    try {
        for(int t_ind = 0; t_ind < num_time_steps; t_ind++) {
            // (1) add the bound of robust input (||v||)
            Interval rho_max_temp = Interval(0.0);
            for (int i = 0; i < num_motors; i++) {
                // compute norm of disturbance
                const PZSparse model_uncertainty_disturbance = 
                    data_sparses_interval[t_ind].tau(i) - data_sparses[t_ind].tau(i);
                const Interval temp = model_uncertainty_disturbance.toInterval();
                rho_max_temp += temp * temp;

                torque_radii(i, t_ind) = 
                    robust_input_bound + 
                    0.5 * std::max(std::abs(temp.lower()), std::abs(temp.upper()));
            }
            rho_max_temp = sqrt(rho_max_temp);
            
            for (int i = 0; i < num_motors; i++) {
                torque_radii(i, t_ind) += 0.5 * rho_max_temp.upper();
            }

            // (2) add friction
            for (int i = 0; i < num_motors; i++) {
                torque_radii(i, t_ind) += robotInfoPtr_->model.friction(i);
            }

            // (3) add the radius back to the nominal input PZ (after reducing)
            for (int i = 0; i < num_motors; i++) {
                data_sparses[t_ind].tau(i).independent += torque_radii(i, t_ind);
            }
            // so that dynPtr_->torque_nom would now store the torque PZs with robust input bounds
        
            // compute the maximum of uncertainty in the sphere centers over all time intervals
            for (int i = 0; i < robotInfoPtr_->num_spheres; i++) {
                const std::string sphere_name = "collision-" + std::to_string(i);

                // if (robotInfoPtr_->model.existFrameName(sphere_name)) {
                    const pinocchio::FrameIndex frame_id = 
                        robotInfoPtr_->model.getFrameId(sphere_name);

                    const auto& sphere_center = data_sparses[t_ind].oMf[frame_id].translation();

                    const double x_uncertainty = sphere_center(0).independent;
                    const double y_uncertainty = sphere_center(1).independent;
                    const double z_uncertainty = sphere_center(2).independent;
                    const double total_uncertainty = std::sqrt(x_uncertainty * x_uncertainty + 
                                                               y_uncertainty * y_uncertainty + 
                                                               z_uncertainty * z_uncertainty);

                    sphere_radii(i, t_ind) = robotInfoPtr_->sphere_radii[i] + total_uncertainty;
                // }
                // else {
                //     throw std::runtime_error("Frame collision-" + std::to_string(num_spheres) + " not found!");
                // }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error computing robust input bounds! " << e.what() << std::endl;
        throw std::runtime_error("Error computing robust input bounds! Check previous error message!");
    }
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR