#include "ReachableSets.h"

namespace RAPTOR {
namespace Armour {

void GenerateJRS(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                 std::shared_ptr<BezierCurveInterval>& trajPtr_) {
    int t_ind = 0; // openmp loop index
    const size_t num_time_steps = trajPtr_->num_time_steps;
    try {
        #pragma omp parallel for shared(trajPtr_) private(t_ind) schedule(dynamic)
        for(t_ind = 0; t_ind < num_time_steps; t_ind++) {
            trajPtr_->makePolyZono(t_ind);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error creating JRS! " << e.what() << std::endl;
        throw std::runtime_error("Error creating JRS! Check previous error message!");
    }               
}

void GenerateLinkAndTorquePZs(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                              std::shared_ptr<BezierCurveInterval>& trajPtr_,
                              std::shared_ptr<KinematicsDynamics>& kdPtr_) {
    int t_ind = 0; // openmp loop index
    const size_t num_time_steps = trajPtr_->num_time_steps;
    try {
        #pragma omp parallel for shared(kdPtr_) private(t_ind) schedule(dynamic)
        for(t_ind = 0; t_ind < num_time_steps; t_ind++) {
            // compute link PZs through forward kinematics
            kdPtr_->fk(t_ind);

            // reduce non-only-k-dependent generators so that slice takes less time
            for (int i = 0; i < 3 * robotInfoPtr_->num_spheres; i++) {
                kdPtr_->sphere_centers(i, t_ind).reduce();
            }

            // compute nominal torque
            kdPtr_->rnea_nominal(t_ind);

            // compute interval torque
            kdPtr_->rnea_interval(t_ind);

            // compute max disturbance (stored in torque_int)
            for (int i = 0; i < NUM_FACTORS; i++) {
                kdPtr_->torque_int(i, t_ind) = kdPtr_->torque_int(i, t_ind) - kdPtr_->torque_nom(i, t_ind);
            }

            // reduce non-only-k-dependent generators so that slice takes less time
            for (int i = 0; i < NUM_FACTORS; i++) {
                kdPtr_->torque_nom(i, t_ind).reduce();
            }

            // reduce non-only-k-dependent generators so that slice takes less time for contact wrench
            for (int i = 0; i < 3 * (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors); i++) {
                kdPtr_->contact_force_int(i, t_ind).reduce();
                kdPtr_->contact_moment_int(i, t_ind).reduce();
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error computing link PZs and nominal torque PZs! " << e.what() << std::endl;
        throw std::runtime_error("Error computing link PZs and nominal torque PZs! Check previous error message!");
    }
}

Eigen::MatrixXd ComputeRobustInputBounds(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                                         std::shared_ptr<BezierCurveInterval>& trajPtr_,
                                         std::shared_ptr<KinematicsDynamics>& kdPtr_) {
    const double robust_input_bound = 
        0.5 * robotInfoPtr_->ultimate_bound_info.alpha * 
            (robotInfoPtr_->ultimate_bound_info.M_max - robotInfoPtr_->ultimate_bound_info.M_min) * 
                robotInfoPtr_->ultimate_bound_info.eps;                            

    const size_t num_time_steps = trajPtr_->num_time_steps;
    const size_t num_motors = robotInfoPtr_->num_motors;

    Eigen::MatrixXd torque_radius(num_motors, num_time_steps);

    try {
        std::vector<double> sphere_center_uncertainty(robotInfoPtr_->num_spheres, 0.0);

        for(int t_ind = 0; t_ind < num_time_steps; t_ind++) {
            // (1) add the bound of robust input (||v||)
            Interval rho_max_temp = Interval(0.0);
            for (int i = 0; i < num_motors; i++) {
                // compute norm of disturbance
                Interval temp = kdPtr_->torque_int(i, t_ind).toInterval(); // should be a 1-dim Interval
                rho_max_temp += temp * temp;

                torque_radius(i, t_ind) = robust_input_bound + 0.5 * std::max(abs(temp.lower()), abs(temp.upper()));
            }
            rho_max_temp = sqrt(rho_max_temp);
            
            for (int i = 0; i < num_motors; i++) {
                torque_radius(i, t_ind) += 0.5 * rho_max_temp.upper();
            }

            // (2) add friction
            for (int i = 0; i < num_motors; i++) {
                torque_radius(i, t_ind) += robotInfoPtr_->model.friction[i];
            }

            // (3) add the radius back to the nominal input PZ (after reducing)
            for (int i = 0; i < num_motors; i++) {
                kdPtr_->torque_nom(i, t_ind).independent += torque_radius(i);
            }
            // so that kdPtr_->torque_nom would now store the torque PZs with robust input bounds`
        
            // compute the maximum of uncertainty in the sphere centers over all time intervals
            for (int i = 0; i < robotInfoPtr_->num_spheres; i++) {
                const double x_uncertainty = kdPtr_->sphere_centers(3 * i + 0, t_ind).independent;
                const double y_uncertainty = kdPtr_->sphere_centers(3 * i + 1, t_ind).independent;
                const double z_uncertainty = kdPtr_->sphere_centers(3 * i + 2, t_ind).independent;
                const double total_uncertainty = sqrtf(x_uncertainty * x_uncertainty + 
                                                      y_uncertainty * y_uncertainty + 
                                                      z_uncertainty * z_uncertainty);

                sphere_center_uncertainty[i] = std::max(sphere_center_uncertainty[i], total_uncertainty);
            }
        }

        // add the maximum uncertainty to the original sphere radii
        for (int i = 0; i < robotInfoPtr_->num_spheres; i++) {
            kdPtr_->sphere_radii[i] += sphere_center_uncertainty[i];
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error computing robust input bounds! " << e.what() << std::endl;
        throw std::runtime_error("Error computing robust input bounds! Check previous error message!");
    }

    return torque_radius;
}

}; // namespace Armour
}; // namespace RAPTOR