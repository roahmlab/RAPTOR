#include "EndEffectorParametersIdentificationMomentum.h"

namespace RAPTOR {

// // constructor
// EndEffectorParametersIdentificationMomentum::EndEffectorParametersIdentificationMomentum()
// {
// }

// // destructors
// EndEffectorParametersIdentificationMomentum::~EndEffectorParametersIdentificationMomentum()
// {
// }

bool EndEffectorParametersIdentificationMomentum::set_parameters(
    const Model& model_input,
    const VecXd offset_input
)
{ 
    enable_hessian = true;

    // parse the robot model
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(*modelPtr_);

    phi = VecXd::Zero(10 * modelPtr_->nv);
    for (Index i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1;
        phi.segment<10>(10 * i) =
            modelPtr_->inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }
    phi_original = phi;

    std::cout << "End effector estimation from URDF file: " << phi_original.tail(10).transpose() << std::endl;

    offset = offset_input;
    if (offset.size() != modelPtr_->nv) { // offset is disabled
        offset = VecXd::Zero(modelPtr_->nv);
    }

    // simply give 0 as initial guess
    x0 = VecXd::Zero(10);

    return true;
}

void EndEffectorParametersIdentificationMomentum::add_trajectory_file(
    const std::string filename_input,
    const SensorNoiseInfo sensor_noise_input,
    const int H_input,
    const TimeFormat time_format,
    const int downsample_rate) {
    // this trajectory is to compute momentum regressors
    trajPtrs_.push_back(std::make_shared<TrajectoryData>(filename_input, 
                                                         sensor_noise_input,
                                                         time_format,
                                                         downsample_rate));

    if (trajPtrs_.back()->Nact != modelPtr_->nv) {
        throw std::invalid_argument("The number of active joints in the trajectory does not match the robot model.");
    }

    trajectoryFilenames_.push_back(filename_input);

    // this trajectory is to compute gravity regressors,
    // so set velocity to 0 while acceleration is already 0 in TrajectoryData
    trajPtrs2_.push_back(std::make_shared<TrajectoryData>(trajPtrs_.back()->T, 
                                                          trajPtrs_.back()->N, 
                                                          trajPtrs_.back()->Nact, 
                                                          false, 
                                                          sensor_noise_input));

    for (Index i = 0; i < trajPtrs_.back()->N; i++) {
        trajPtrs2_.back()->tspan(i) = trajPtrs_.back()->tspan(i);
        trajPtrs2_.back()->q(i) = trajPtrs_.back()->q(i);
    }

    initialize_regressors(trajPtrs_.back(),
                          trajPtrs2_.back(),
                          H_input);

    // combine all the regressors
    num_segments.clear();
    total_num_segments = 0;
    for (Index i = 0; i < Aseg.size(); i++) {
        num_segments.push_back(Aseg[i].rows() / modelPtr_->nv);
        total_num_segments += num_segments.back();
    }
}

void EndEffectorParametersIdentificationMomentum::initialize_regressors(const std::shared_ptr<TrajectoryData>& trajPtr_,
                                                                        const std::shared_ptr<TrajectoryData>& trajPtr2_,
                                                                        const int H_input) {
    // create regressor compute object
    mrPtr_ = std::make_shared<MomentumRegressor>(*modelPtr_, trajPtr_);
    ridPtr_ = std::make_shared<RegressorInverseDynamics>(*modelPtr_, trajPtr2_, false);

    // forward integration horizon
    H = H_input;
    int num_segment = trajPtr_->N / H - 1;

    if (num_segment <= 0) {
        THROW_EXCEPTION(IpoptException, "0 segments");
    }

    // trajectoryData does not require any decision variable, 
    // so simply put an empty Eigen vector as placeholder here
    mrPtr_->compute(VecXd::Zero(1), false);
    ridPtr_->compute(VecXd::Zero(1), false);

    // now compute regression elements A and b
    // which are essentially combination of different dynamic regressors
    MatXd A_seg_i = MatXd::Zero(modelPtr_->nv * num_segment, 10 * modelPtr_->nv);
    VecXd b_seg_i = VecXd::Zero(modelPtr_->nv * num_segment);

    Index i = 0;
    #pragma omp parallel for shared(modelPtr_, trajPtr_,  mrPtr_, ridPtr_, A_seg_i, b_seg_i) private(i) schedule(dynamic, 1)
    for (i = 0; i < num_segment; i++) {
        const int seg_start = i * H;
        const int seg_end = seg_start + H;

        const MatXd& Y_Hqd_1 = mrPtr_->Y.middleRows(seg_start * modelPtr_->nv, modelPtr_->nv);
        const MatXd& Y_Hqd_2 = mrPtr_->Y.middleRows(seg_end * modelPtr_->nv, modelPtr_->nv);

        MatXd int_Y_CTqd_g = MatXd::Zero(modelPtr_->nv, 10 * modelPtr_->nv);
        VecXd int_ctrl = VecXd::Zero(modelPtr_->nv);

        for (int j = seg_start; j < seg_end; j++) {
            const double dt = trajPtr_->tspan(j+1) - trajPtr_->tspan(j);

            const MatXd& Y_CTv_i = mrPtr_->Y_CTv.middleRows(j * modelPtr_->nv, modelPtr_->nv);
            const MatXd& Yg_i = ridPtr_->Y.middleRows(j * modelPtr_->nv, modelPtr_->nv);

            int_Y_CTqd_g += (Y_CTv_i - Yg_i) * dt;

            // Note that here trajPtr_->q_dd stores the applied torque
            int_ctrl += (trajPtr_->q_dd(j) -
                         modelPtr_->friction.cwiseProduct(trajPtr_->q_d(j).cwiseSign()) -
                         modelPtr_->damping.cwiseProduct(trajPtr_->q_d(j)) -
                         offset) * dt;
        }

        A_seg_i.middleRows(i * modelPtr_->nv, modelPtr_->nv) = (Y_Hqd_2 - Y_Hqd_1) - int_Y_CTqd_g;
        b_seg_i.segment(i * modelPtr_->nv, modelPtr_->nv) = int_ctrl - 
                                                            modelPtr_->armature.cwiseProduct(
                                                                trajPtr_->q_d(seg_end) - trajPtr_->q_d(seg_start));
    }

    Aseg.push_back(A_seg_i);
    bseg.push_back(b_seg_i);
}

void EndEffectorParametersIdentificationMomentum::reset() {
    Optimizer::reset();

    trajPtrs_.clear();
    trajPtrs2_.clear();
    trajectoryFilenames_.clear();

    Aseg.clear();
    bseg.clear();

    num_segments.clear();

    std::cout << "End effector identification reset, all trajectory data removed." << std::endl;
}

void EndEffectorParametersIdentificationMomentum::finalize_solution(
    SolverReturn               status,
    Index                      n,
    const Number*              x,
    const Number*              z_L,
    const Number*              z_U,
    Index                      m,
    const Number*              g,
    const Number*              lambda,
    Number                     obj_value,
    const IpoptData*           ip_data,
    IpoptCalculatedQuantities* ip_cq
) {
    Optimizer::finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);

    theta_solution = z_to_theta(solution);

    std::cout << "Performing error analysis" << std::endl;

    MatXd A(modelPtr_->nv * total_num_segments, 10 * modelPtr_->nv);
    VecXd b(modelPtr_->nv * total_num_segments);
    A.setZero();
    b.setZero();
    Index row_start = 0;
    for (Index i = 0; i < Aseg.size(); i++) {
        const MatX& Aseg_i = Aseg[i];
        const MatX& bseg_i = bseg[i];
        A.middleRows(row_start, Aseg_i.rows()) = Aseg_i;
        b.segment(row_start, bseg_i.size()) = bseg_i;
        row_start += Aseg_i.rows();
    }

    const MatXd& A_opt = A.rightCols(10);
    const VecXd b_opt = b - A * phi_original;

    Mat10d dtheta;
    Eigen::Array<Mat10d, 1, 10> ddtheta;
    phi.tail(10) = dd_z_to_theta(solution, dtheta, ddtheta);
    VecXd diff = A * phi - b;

    MatXd temp1 = A_opt * dtheta;
    MatXd temp2 = diff.transpose() * A_opt;
    Mat10d temp3 = temp1.transpose() * temp1;
    Mat10d temp4;
    temp4.setZero();
    for (Index i = 0; i < 10; i++) {
        temp4 += temp2(i) * ddtheta(i);
    }
    Mat10d p_z_p_eta = temp3 + temp4;
    Eigen::LDLT<MatXd> ldlt(p_z_p_eta);
    Mat10d p_z_p_eta_inv;
    if (ldlt.info() == Eigen::Success) {
        p_z_p_eta_inv = ldlt.solve(Mat10d::Identity());
        // std::cout << dtheta * p_z_p_eta_inv << std::endl;
    }
    else {
        std::cerr << p_z_p_eta << std::endl;
        THROW_EXCEPTION(IpoptException, "*** Chelosky decomposition not successful!");
    }

    theta_uncertainty.setZero();
    Vec10d p_z_p_x;
    Vec10d p_eta_p_x;
    Vec10d p_theta_p_x;

    // compute p_b_p_x for each of the x
    // (1) applied torque data: num_segment * H * modelPtr_->nv
    Index pivot = 0;
    for (Index tid = 0; tid < trajPtrs_.size(); tid++) {
        const int num_segment = num_segments[tid];
        const auto& trajPtr_ = trajPtrs_[tid];
        const auto& sensor_noise = trajPtr_->sensor_noise;
        for (Index s = 0; s < num_segment; s++) {
            int seg_start = s * H;
            int seg_end = seg_start + H;
            for (int j = seg_start; j < seg_end; j++) {
                double dt = trajPtr_->tspan(j + 1) - trajPtr_->tspan(j);
                for (Index k = 0; k < modelPtr_->nv; k++) {
                    // p_b_p_x = dt;
                    p_z_p_x = -dt * temp1.row((pivot + s) * modelPtr_->nv + k);
                    p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                    p_theta_p_x = dtheta * p_eta_p_x;
                    double torque_error = 0.0;
                    if (sensor_noise.acceleration_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
                        torque_error = std::abs(trajPtr_->q_dd(j)(k) * sensor_noise.acceleration_error(k));
                    }
                    else {
                        torque_error = sensor_noise.acceleration_error(k);
                    }
                    theta_uncertainty += p_theta_p_x.cwiseAbs() * torque_error;
                }
            }
        }
        pivot += num_segment;
    }

    // (2) friction parameters: 4 * modelPtr_->nv
    pivot = 0;
    for (Index tid = 0; tid < trajPtrs_.size(); tid++) {
        const int num_segment = num_segments[tid];
        const auto& trajPtr_ = trajPtrs_[tid];
        for (Index s = 0; s < num_segment; s++) {
            int seg_start = s * H;
            int seg_end = seg_start + H;
            for (int j = seg_start; j < seg_end; j++) {
                double dt = trajPtr_->tspan(j + 1) - trajPtr_->tspan(j);
                for (Index k = 0; k < modelPtr_->nv; k++) {
                    // friction
                    // p_b_p_x = dt * Utils::sign(trajPtr_->q_d(j)(k));
                    p_z_p_x = -dt * Utils::sign(trajPtr_->q_d(j)(k)) * temp1.row((pivot + s) * modelPtr_->nv + k);
                    p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                    p_theta_p_x = dtheta * p_eta_p_x;
                    double friction_parameter_error = std::abs(0.05 * modelPtr_->friction(k));
                    theta_uncertainty += p_theta_p_x.cwiseAbs() * friction_parameter_error;

                    // damping
                    // p_b_p_x = dt * trajPtr_->q_d(j)(k);
                    p_z_p_x = -dt * trajPtr_->q_d(j)(k) * temp1.row((pivot + s) * modelPtr_->nv + k);
                    p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                    p_theta_p_x = dtheta * p_eta_p_x;
                    double damping_parameter_error = std::abs(0.05 * modelPtr_->damping(k));
                    theta_uncertainty += p_theta_p_x.cwiseAbs() * damping_parameter_error;

                    // offset
                    // p_b_p_x = dt;
                    p_z_p_x = -dt * temp1.row((pivot + s) * modelPtr_->nv + k);
                    p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                    p_theta_p_x = dtheta * p_eta_p_x;
                    double offset_error = std::abs(0.05 * offset(k));
                    theta_uncertainty += p_theta_p_x.cwiseAbs() * offset_error;
                }
            }
        }
        pivot += num_segment;
    }

    // (3) other link inertial parameters: 10 * (modelPtr_->nv - 1)
    MatXd p_b_p_x = -A.leftCols(10 * (modelPtr_->nv - 1));
    MatXd p_z_p_x_2 = -p_b_p_x.transpose() * temp1;
    MatXd p_eta_p_x_2 = -p_z_p_eta_inv * p_z_p_x_2.transpose();
    MatXd p_theta_p_x_2 = dtheta * p_eta_p_x_2;
    theta_uncertainty += p_theta_p_x_2.cwiseAbs() * phi_original.head(10 * (modelPtr_->nv - 1)).cwiseAbs() * 0.05;

    std::cout << "Uncertainty on the estimated end-effector inertial parameters: " << theta_uncertainty.transpose() << std::endl;
}

}; // namespace RAPTOR
