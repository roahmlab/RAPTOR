#include "EndEffectorParametersIdentification.h"

namespace RAPTOR {

// // constructor
// EndEffectorParametersIdentification::EndEffectorParametersIdentification()
// {
// }

// // destructors
// EndEffectorParametersIdentification::~EndEffectorParametersIdentification()
// {
// }

bool EndEffectorParametersIdentification::set_parameters(
    const Model& model_input,
    const std::string filename_input,
    const SensorNoiseInfo sensor_noise_input,
    const int H_input,
    const TimeFormat time_format,
    const int downsample_rate,
    const VecXd offset_input
)
{ 
    // macro NUM_THREADS should be define in cmake
    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        THROW_EXCEPTION(IpoptException, "macro NUM_THREADS is not defined!");
    #endif

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

    offset = offset_input;
    if (offset.size() != modelPtr_->nv) { // offset is disabled
        offset = VecXd::Zero(modelPtr_->nv);
    }

    // this trajectory is to compute momentum regressors
    trajPtr_ = std::make_shared<TrajectoryData>(filename_input, 
                                                sensor_noise_input,
                                                time_format,
                                                downsample_rate);

    // this trajectory is to compute gravity regressors,
    // so set velocity to 0 while acceleration is already 0 in TrajectoryData
    trajPtr2_ = std::make_shared<TrajectoryData>(filename_input, 
                                                 sensor_noise_input,
                                                 time_format,
                                                 downsample_rate);
    for (Index i = 0; i < trajPtr2_->N; i++) {
        trajPtr2_->q_d(i).setZero();
        trajPtr2_->q_dd(i).setZero();
    }

    // create regressor compute object
    mrPtr_ = std::make_shared<MomentumRegressor>(*modelPtr_, trajPtr_);
    ridPtr_ = std::make_shared<RegressorInverseDynamics>(*modelPtr_, trajPtr2_, false);

    // forward integration horizon
    H = H_input;
    num_segment = trajPtr_->N / H;

    if (num_segment <= 0) {
        THROW_EXCEPTION(IpoptException, "0 segments");
    }

    // trajectoryData does not require any decision variable, 
    // so simply put an empty Eigen vector as placeholder here
    mrPtr_->compute(VecXd::Zero(1), false);
    ridPtr_->compute(VecXd::Zero(1), false);

    // now compute regression elements A and b
    // which are essentially combination of different dynamic regressors
    A.resize(modelPtr_->nv * num_segment, 10 * modelPtr_->nv);
    b.resize(modelPtr_->nv * num_segment);

    Index i = 0;
    #pragma omp parallel for shared(modelPtr_, trajPtr_,  mrPtr_, ridPtr_, A, b) private(i) schedule(dynamic, 1)
    for (i = 0; i < trajPtr_->N - H - 1; i += H) {
        int seg_start = i;
        int seg_end = seg_start + H;

        const MatXd& Y_Hqd_1 = mrPtr_->Y.middleRows(seg_start * modelPtr_->nv, modelPtr_->nv);
        const MatXd& Y_Hqd_2 = mrPtr_->Y.middleRows(seg_end * modelPtr_->nv, modelPtr_->nv);

        MatXd int_Y_CTqd_g = MatXd::Zero(modelPtr_->nv, 10 * modelPtr_->nv);
        VecXd int_ctrl = VecXd::Zero(modelPtr_->nv);

        for (int j = seg_start; j < seg_end; j++) {
            double dt = trajPtr_->tspan(j+1) - trajPtr_->tspan(j);

            MatXd Y_CTqd_i = mrPtr_->Y_CTv.middleRows(j * modelPtr_->nv, modelPtr_->nv);
            MatXd Yg_i = ridPtr_->Y.middleRows(j * modelPtr_->nv, modelPtr_->nv);

            int_Y_CTqd_g += (Y_CTqd_i - Yg_i) * dt;

            // Note that here trajPtr_->q_dd stores the applied torque
            int_ctrl += (trajPtr_->q_dd(j) -
                         modelPtr_->friction.cwiseProduct(trajPtr_->q_d(j).cwiseSign()) -
                         modelPtr_->damping.cwiseProduct(trajPtr_->q_d(j)) -
                         offset) * dt;
        }
        
        int s = i / H;
        A.middleRows(s * modelPtr_->nv, modelPtr_->nv) = (Y_Hqd_2 - Y_Hqd_1) - int_Y_CTqd_g;
        b.segment(s * modelPtr_->nv, modelPtr_->nv) = int_ctrl - modelPtr_->armature.cwiseProduct(trajPtr_->q_d(seg_end) - trajPtr_->q_d(seg_start));
    }

    // simply give 0 as initial guess
    x0 = VecXd::Zero(10); 

    // parse sensor noise
    if (sensor_noise_input.position_error.norm() == 0.0 &&
        sensor_noise_input.velocity_error.norm() == 0.0 &&
        sensor_noise_input.acceleration_error.norm() == 0.0) {
        std::cout << "No sensor noise is added" << std::endl;
    }
    else {
        std::cout << "Sensor noise is added" << std::endl;
        
        // mrIntPtr_ = std::make_shared<IntervalMomentumRegressor>(*modelPtr_, trajPtr_);
        // ridIntPtr_ = std::make_shared<IntervalRegressorInverseDynamics>(*modelPtr_, trajPtr2_);
    }

    return true;
}

bool EndEffectorParametersIdentification::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n = 10;       // End-effector parameters 
    numVars = n;

    // number of inequality constraint
    numCons = 0;
    m = numCons;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool EndEffectorParametersIdentification::eval_f(
    Index n,
    const Number *x,
    bool new_x,
    Number &obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecXd z = Utils::initializeEigenVectorFromArray(x, n);

    // Update the inertial parameters
    phi.tail(10) = z_to_theta(z);

    // Compute the ojective function
    VecXd diff = A * phi - b;

    obj_value = 0.5 * diff.dot(diff) / b.size();

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool EndEffectorParametersIdentification::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecXd z = Utils::initializeEigenVectorFromArray(x, n);
    VecXd grad_f_vec = VecXd::Zero(n);
    
    Mat10d dtheta;
    phi.tail(10) = d_z_to_theta(z, dtheta);
    VecXd diff = A * phi - b;

    // Compute the gradient
    grad_f_vec = (diff.transpose() * A.rightCols(10) * dtheta) / b.size();

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i);
    }

    return true;
}

bool EndEffectorParametersIdentification::eval_hess_f(
    Index         n,
    const Number* x,
    bool          new_x,
    MatX&         hess_f
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    VecXd z = Utils::initializeEigenVectorFromArray(x, n);
    hess_f = MatX::Zero(n, n);

    Mat10d dtheta;
    Eigen::Array<Mat10d, 1, 10> ddtheta;
    phi.tail(10) = dd_z_to_theta(z, dtheta, ddtheta);
    VecXd diff = A * phi - b;

    // Compute the Hessian
    MatX temp1 = A.rightCols(10) * dtheta;
    hess_f = temp1.transpose() * temp1;

    MatX temp2 = diff.transpose() * A.rightCols(10);
    for (Index i = 0; i < n; i++) {
        hess_f += temp2(i) * ddtheta(i);
    }

    hess_f /= b.size();

    return true;
}

void EndEffectorParametersIdentification::finalize_solution(
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

    const auto& sensor_noise = trajPtr_->sensor_noise;

    Mat10d dtheta;
    Eigen::Array<Mat10d, 1, 10> ddtheta;
    phi.tail(10) = dd_z_to_theta(solution, dtheta, ddtheta);
    VecXd diff = A * phi - b;

    MatXd temp1 = A.rightCols(10) * dtheta;
    MatXd temp2 = diff.transpose() * A;
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
        std::cout << p_z_p_eta_inv << std::endl;
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
    for (Index i = 0; i < trajPtr_->N - H - 1; i += H) {
        int seg_start = i;
        int seg_end = seg_start + H;
        int s = i / H;
        for (int j = seg_start; j < seg_end; j++) {
            double dt = trajPtr_->tspan(j + 1) - trajPtr_->tspan(j);
            for (Index k = 0; k < modelPtr_->nv; k++) {
                // p_b_p_x = dt;
                p_z_p_x = -dt * temp1.row(s * modelPtr_->nv + k);
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

    // (2) friction parameters: 4 * modelPtr_->nv
    for (Index i = 0; i < trajPtr_->N - H - 1; i += H) {
        int seg_start = i;
        int seg_end = seg_start + H;
        int s = i / H;
        for (int j = seg_start; j < seg_end; j++) {
            double dt = trajPtr_->tspan(j + 1) - trajPtr_->tspan(j);
            for (Index k = 0; k < modelPtr_->nv; k++) {
                // friction
                // p_b_p_x = dt * Utils::sign(trajPtr_->q_d(j)(k));
                p_z_p_x = -dt * Utils::sign(trajPtr_->q_d(j)(k)) * temp1.row(s * modelPtr_->nv + k);
                p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                p_theta_p_x = dtheta * p_eta_p_x;
                double friction_parameter_error = std::abs(0.05 * modelPtr_->friction(k));
                theta_uncertainty += p_theta_p_x.cwiseAbs() * friction_parameter_error;

                // damping
                // p_b_p_x = dt * trajPtr_->q_d(j)(k);
                p_z_p_x = -dt * trajPtr_->q_d(j)(k) * temp1.row(s * modelPtr_->nv + k);
                p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                p_theta_p_x = dtheta * p_eta_p_x;
                double damping_parameter_error = std::abs(0.05 * modelPtr_->damping(k));
                theta_uncertainty += p_theta_p_x.cwiseAbs() * damping_parameter_error;

                // offset
                // p_b_p_x = dt;
                p_z_p_x = -dt * temp1.row(s * modelPtr_->nv + k);
                p_eta_p_x = -p_z_p_eta_inv * p_z_p_x;
                p_theta_p_x = dtheta * p_eta_p_x;
                double offset_error = std::abs(0.05 * offset(k));
                theta_uncertainty += p_theta_p_x.cwiseAbs() * offset_error;
            }
        }
    }

    // (3) other link inertial parameters: 10 * (modelPtr_->nv - 1)
    MatXd p_b_p_x = -A.leftCols(10 * (modelPtr_->nv - 1));
    MatXd p_z_p_x_2 = -p_b_p_x.transpose() * temp1;
    MatXd p_eta_p_x_2 = -p_z_p_eta_inv * p_z_p_x_2.transpose();
    MatXd p_theta_p_x_2 = dtheta * p_eta_p_x_2;
    theta_uncertainty += p_theta_p_x_2.cwiseAbs() * phi_original.head(10 * (modelPtr_->nv - 1)).cwiseAbs() * 0.05;

    std::cout << "Uncertainty on the estimated end-effector inertial parameters: " << theta_uncertainty.transpose() << std::endl;
}

Eigen::Vector<double, 10> EndEffectorParametersIdentification::z_to_theta(const VecXd& z) {
    // log-Cholesky parameterization
    // const double d1 = z[0];
    // const double d2 = z[1];
    // const double d3 = z[2];
    // const double d4 = z[3];
    // const double s12 = z[4];
    // const double s23 = z[5];
    // const double s13 = z[6];
    // const double t1 = z[7];
    // const double t2 = z[8];
    // const double t3 = z[9];

    // // This is the direct way to compute theta by definition of Log Cheloysky decomposition
    // Mat4d U;
    // U << std::exp(d1), s12,          s13,          t1,
    //      0.0,          std::exp(d2), s23,          t2,
    //      0.0,          0.0,          std::exp(d3), t3,
    //      0.0,          0.0,          0.0,          std::exp(d4);

    // // Compute LMI = U' * U
    // Mat4d LMI = U.transpose() * U;

    // // End-effector z
    // VecXd theta = VecXd::Zero(10);
    // theta(0) = LMI(3, 3);                        
    // theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    // theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    // theta(5) = -LMI(0, 1);                       // IXY
    // theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    // theta(7) = -LMI(0, 2);                       // IXZ
    // theta(8) = -LMI(1, 2);                       // IYZ
    // theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    Vec10d theta;
    // double t5 = std::exp(d1);
    // double t6 = std::exp(d2);
    // double t7 = d1*2.0;
    // double t8 = d2*2.0;
    // double t9 = d3*2.0;
    // double t10 = s12*s12;
    // double t11 = s13*s13;
    // double t12 = s23*s23;
    // double t13 = std::exp(t7);
    // double t14 = std::exp(t8);
    // double t15 = std::exp(t9);
    // theta(0) = std::exp(d4*2.0)+t1*t1+t2*t2+t3*t3;
    // theta(1) = t1*t5;
    // theta(2) = s12*t1+t2*t6;
    // theta(3) = s13*t1+s23*t2+t3*std::exp(d3);
    // theta(4) = t10+t11+t12+t14+t15;
    // theta(5) = -s12*t5;
    // theta(6) = t11+t12+t13+t15;
    // theta(7) = -s13*t5;
    // theta(8) = -s12*s13-s23*t6;
    // theta(9) = t10+t13+t14;
    
    const double alpha = z[0];
    const double d1 = z[1];
    const double d2 = z[2];
    const double d3 = z[3];
    const double s12 = z[4];
    const double s23 = z[5];
    const double s13 = z[6];
    const double t1 = z[7];
    const double t2 = z[8];
    const double t3 = z[9];
    const double exp_d1 = std::exp(d1);
    const double exp_d2 = std::exp(d2);
    const double exp_d3 = std::exp(d3);
    theta[0] = 1;
    theta[1] = t1;
    theta[2] = t2;
    theta[3] = t3;
    theta[4] = std::pow(s23, 2) + std::pow(t2, 2) + std::pow(t3, 2) + std::pow(exp_d2, 2) + std::pow(exp_d3, 2);
    theta[5] = -s12 * exp_d2 - s13 * s23 - t1 * t2;
    theta[6] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(t1, 2) + std::pow(t3, 2) + std::pow(exp_d1, 2) + std::pow(exp_d3, 2);
    theta[7] = -s13 * exp_d3 - t1 * t3;
    theta[8] = -s23 * exp_d3 - t2 * t3;
    theta[9] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(s23, 2) + std::pow(t1, 2) + std::pow(t2, 2) + std::pow(exp_d1, 2) + std::pow(exp_d2, 2);
    const double exp_2_alpha = std::exp(2 * alpha);
    theta *= exp_2_alpha;

    return theta;
}

Eigen::Vector<double, 10> EndEffectorParametersIdentification::d_z_to_theta(
    const VecXd& z,
    Mat10d& dtheta) {
    // log-Cholesky parameterization
    // const double d1 = z[0];
    // const double d2 = z[1];
    // const double d3 = z[2];
    // const double d4 = z[3];
    // const double s12 = z[4];
    // const double s23 = z[5];
    // const double s13 = z[6];
    // const double t1 = z[7];
    // const double t2 = z[8];
    // const double t3 = z[9];

    // Vec10d theta;
    // double t5 = std::exp(d1);
    // double t6 = std::exp(d2);
    // double t7 = d1*2.0;
    // double t8 = d2*2.0;
    // double t9 = d3*2.0;
    // double t10 = s12*s12;
    // double t11 = s13*s13;
    // double t12 = s23*s23;
    // double t13 = std::exp(t7);
    // double t14 = std::exp(t8);
    // double t15 = std::exp(t9);
    // theta(0) = std::exp(d4*2.0)+t1*t1+t2*t2+t3*t3;
    // theta(1) = t1*t5;
    // theta(2) = s12*t1+t2*t6;
    // theta(3) = s13*t1+s23*t2+t3*std::exp(d3);
    // theta(4) = t10+t11+t12+t14+t15;
    // theta(5) = -s12*t5;
    // theta(6) = t11+t12+t13+t15;
    // theta(7) = -s13*t5;
    // theta(8) = -s12*s13-s23*t6;
    // theta(9) = t10+t13+t14;

    // t7 = std::exp(d3);
    // t8 = d1*2.0;
    // t9 = d2*2.0;
    // t10 = d3*2.0;
    // t11 = s12*2.0;
    // t12 = s13*2.0;
    // t13 = s23*2.0;
    // t14 = std::exp(t8);
    // t15 = std::exp(t9);
    // double t16 = std::exp(t10);
    // double t17 = -t5;
    // double t18 = t14*2.0;
    // double t19 = t15*2.0;
    // double t20 = t16*2.0;

    // dtheta.setZero();
    // dtheta(0, 3) = std::exp(d4 * 2.0) * 2.0;
    // dtheta(0, 7) = t1 * 2.0;
    // dtheta(0, 8) = t2 * 2.0;
    // dtheta(0, 9) = t3 * 2.0;
    // dtheta(1, 0) = t1 * t5;
    // dtheta(1, 7) = t5;
    // dtheta(2, 1) = t2 * t6;
    // dtheta(2, 4) = t1;
    // dtheta(2, 7) = s12;
    // dtheta(2, 8) = t6;
    // dtheta(3, 2) = t3 * t7;
    // dtheta(3, 5) = t2;
    // dtheta(3, 6) = t1;
    // dtheta(3, 7) = s13;
    // dtheta(3, 8) = s23;
    // dtheta(3, 9) = t7;
    // dtheta(4, 1) = t19;
    // dtheta(4, 2) = t20;
    // dtheta(4, 4) = t11;
    // dtheta(4, 5) = t13;
    // dtheta(4, 6) = t12;
    // dtheta(5, 0) = s12 * t17;
    // dtheta(5, 4) = t17;
    // dtheta(6, 0) = t18;
    // dtheta(6, 2) = t20;
    // dtheta(6, 5) = t13;
    // dtheta(6, 6) = t12;
    // dtheta(7, 0) = s13 * t17;
    // dtheta(7, 6) = t17;
    // dtheta(8, 1) = -s23 * t6;
    // dtheta(8, 4) = -s13;
    // dtheta(8, 5) = -t6;
    // dtheta(8, 6) = -s12;
    // dtheta(9, 0) = t18;
    // dtheta(9, 1) = t19;
    // dtheta(9, 4) = t11;

    Vec10d theta;

    const double alpha = z[0];
    const double d1 = z[1];
    const double d2 = z[2];
    const double d3 = z[3];
    const double s12 = z[4];
    const double s23 = z[5];
    const double s13 = z[6];
    const double t1 = z[7];
    const double t2 = z[8];
    const double t3 = z[9];

    const double exp_2alpha = std::exp(2 * alpha);
    const double exp_2d1 = std::exp(2 * d1);
    const double exp_2d2 = std::exp(2 * d2);
    const double exp_2d3 = std::exp(2 * d3);
    const double exp_d1 = std::exp(d1);
    const double exp_d2 = std::exp(d2);
    const double exp_d3 = std::exp(d3);

    theta[0] = 1;
    theta[1] = t1;
    theta[2] = t2;
    theta[3] = t3;
    theta[4] = std::pow(s23, 2) + std::pow(t2, 2) + std::pow(t3, 2) + std::pow(exp_d2, 2) + std::pow(exp_d3, 2);
    theta[5] = -s12 * exp_d2 - s13 * s23 - t1 * t2;
    theta[6] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(t1, 2) + std::pow(t3, 2) + std::pow(exp_d1, 2) + std::pow(exp_d3, 2);
    theta[7] = -s13 * exp_d3 - t1 * t3;
    theta[8] = -s23 * exp_d3 - t2 * t3;
    theta[9] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(s23, 2) + std::pow(t1, 2) + std::pow(t2, 2) + std::pow(exp_d1, 2) + std::pow(exp_d2, 2);
    const double exp_2_alpha = std::exp(2 * alpha);
    theta *= exp_2_alpha;

    dtheta.setZero();

    dtheta(0, 0) = 2 * exp_2alpha;

    dtheta(1, 0) = 2 * t1 * exp_2alpha;
    dtheta(1, 7) = exp_2alpha;

    dtheta(2, 0) = 2 * t2 * exp_2alpha;
    dtheta(2, 8) = exp_2alpha;

    dtheta(3, 0) = 2 * t3 * exp_2alpha;
    dtheta(3, 9) = exp_2alpha;

    dtheta(4, 0) = 2 * (std::pow(s23, 2) + std::pow(t2, 2) + std::pow(t3, 2) + exp_2d2 + exp_2d3) * exp_2alpha;
    dtheta(4, 2) = 2 * exp_2alpha * exp_2d2;
    dtheta(4, 3) = 2 * exp_2alpha * exp_2d3;
    dtheta(4, 5) = 2 * s23 * exp_2alpha;
    dtheta(4, 8) = 2 * t2 * exp_2alpha;
    dtheta(4, 9) = 2 * t3 * exp_2alpha;

    dtheta(5, 0) = -2 * (s12 * exp_d2 + s13 * s23 + t1 * t2) * exp_2alpha;
    dtheta(5, 2) = -s12 * exp_2alpha * exp_d2;
    dtheta(5, 4) = -exp_2alpha * exp_d2;
    dtheta(5, 5) = -s13 * exp_2alpha;
    dtheta(5, 6) = -s23 * exp_2alpha;
    dtheta(5, 7) = -t2 * exp_2alpha;
    dtheta(5, 8) = -t1 * exp_2alpha;

    dtheta(6, 0) = 2 * (std::pow(s12, 2) + std::pow(s13, 2) + std::pow(t1, 2) + std::pow(t3, 2) + exp_2d1 + exp_2d3) * exp_2alpha;
    dtheta(6, 1) = 2 * exp_2alpha * exp_2d1;
    dtheta(6, 3) = 2 * exp_2alpha * exp_2d3;
    dtheta(6, 4) = 2 * s12 * exp_2alpha;
    dtheta(6, 6) = 2 * s13 * exp_2alpha;
    dtheta(6, 7) = 2 * t1 * exp_2alpha;
    dtheta(6, 9) = 2 * t3 * exp_2alpha;

    dtheta(7, 0) = -2 * (s13 * exp_d3 + t1 * t3) * exp_2alpha;
    dtheta(7, 3) = -s13 * exp_2alpha * exp_d3;
    dtheta(7, 6) = -exp_2alpha * exp_d3;
    dtheta(7, 7) = -t3 * exp_2alpha;
    dtheta(7, 9) = -t1 * exp_2alpha;

    dtheta(8, 0) = -2 * (s23 * exp_d3 + t2 * t3) * exp_2alpha;
    dtheta(8, 3) = -s23 * exp_2alpha * exp_d3;
    dtheta(8, 5) = -exp_2alpha * exp_d3;
    dtheta(8, 8) = -t3 * exp_2alpha;
    dtheta(8, 9) = -t2 * exp_2alpha;

    dtheta(9, 0) = 2 * (std::pow(s12, 2) + std::pow(s13, 2) + std::pow(s23, 2) + std::pow(t1, 2) + std::pow(t2, 2) + exp_2d1 + exp_2d2) * exp_2alpha;
    dtheta(9, 1) = 2 * exp_2alpha * exp_2d1;
    dtheta(9, 2) = 2 * exp_2alpha * exp_2d2;
    dtheta(9, 4) = 2 * s12 * exp_2alpha;
    dtheta(9, 5) = 2 * s23 * exp_2alpha;
    dtheta(9, 6) = 2 * s13 * exp_2alpha;
    dtheta(9, 7) = 2 * t1 * exp_2alpha;
    dtheta(9, 8) = 2 * t2 * exp_2alpha;

    return theta;
}

Eigen::Vector<double, 10> EndEffectorParametersIdentification::dd_z_to_theta(
    const VecXd& z,
    Mat10d& dtheta,
    Eigen::Array<Mat10d, 1, 10>& ddtheta) {
    // // log-Cholesky parameterization
    // const double d1 = z[0];
    // const double d2 = z[1];
    // const double d3 = z[2];
    // const double d4 = z[3];
    // const double s12 = z[4];
    // const double s23 = z[5];
    // const double s13 = z[6];
    // const double t1 = z[7];
    // const double t2 = z[8];
    // const double t3 = z[9];

    // Vec10d theta;
    // double t5 = std::exp(d1);
    // double t6 = std::exp(d2);
    // double t7 = d1*2.0;
    // double t8 = d2*2.0;
    // double t9 = d3*2.0;
    // double t10 = s12*s12;
    // double t11 = s13*s13;
    // double t12 = s23*s23;
    // double t13 = std::exp(t7);
    // double t14 = std::exp(t8);
    // double t15 = std::exp(t9);
    // theta(0) = std::exp(d4*2.0)+t1*t1+t2*t2+t3*t3;
    // theta(1) = t1*t5;
    // theta(2) = s12*t1+t2*t6;
    // theta(3) = s13*t1+s23*t2+t3*std::exp(d3);
    // theta(4) = t10+t11+t12+t14+t15;
    // theta(5) = -s12*t5;
    // theta(6) = t11+t12+t13+t15;
    // theta(7) = -s13*t5;
    // theta(8) = -s12*s13-s23*t6;
    // theta(9) = t10+t13+t14;

    // t7 = std::exp(d3);
    // t8 = d1*2.0;
    // t9 = d2*2.0;
    // t10 = d3*2.0;
    // t11 = s12*2.0;
    // t12 = s13*2.0;
    // t13 = s23*2.0;
    // t14 = std::exp(t8);
    // t15 = std::exp(t9);
    // double t16 = std::exp(t10);
    // double t17 = -t5;
    // double t18 = t14*2.0;
    // double t19 = t15*2.0;
    // double t20 = t16*2.0;

    // dtheta.setZero();
    // dtheta(0, 3) = std::exp(d4 * 2.0) * 2.0;
    // dtheta(0, 7) = t1 * 2.0;
    // dtheta(0, 8) = t2 * 2.0;
    // dtheta(0, 9) = t3 * 2.0;
    // dtheta(1, 0) = t1 * t5;
    // dtheta(1, 7) = t5;
    // dtheta(2, 1) = t2 * t6;
    // dtheta(2, 4) = t1;
    // dtheta(2, 7) = s12;
    // dtheta(2, 8) = t6;
    // dtheta(3, 2) = t3 * t7;
    // dtheta(3, 5) = t2;
    // dtheta(3, 6) = t1;
    // dtheta(3, 7) = s13;
    // dtheta(3, 8) = s23;
    // dtheta(3, 9) = t7;
    // dtheta(4, 1) = t19;
    // dtheta(4, 2) = t20;
    // dtheta(4, 4) = t11;
    // dtheta(4, 5) = t13;
    // dtheta(4, 6) = t12;
    // dtheta(5, 0) = s12 * t17;
    // dtheta(5, 4) = t17;
    // dtheta(6, 0) = t18;
    // dtheta(6, 2) = t20;
    // dtheta(6, 5) = t13;
    // dtheta(6, 6) = t12;
    // dtheta(7, 0) = s13 * t17;
    // dtheta(7, 6) = t17;
    // dtheta(8, 1) = -s23 * t6;
    // dtheta(8, 4) = -s13;
    // dtheta(8, 5) = -t6;
    // dtheta(8, 6) = -s12;
    // dtheta(9, 0) = t18;
    // dtheta(9, 1) = t19;
    // dtheta(9, 4) = t11;

    // for (Index i = 0; i < 10; i++) {
    //     ddtheta(i).setZero();
    // }

    // t11 = exp(t8);
    // t12 = exp(t9);
    // t13 = exp(t10);
    // t14 = -t5;
    // t15 = -t6;
    // t16 = t11*4.0;
    // t17 = t12*4.0;
    // t18 = t13*4.0;
    // ddtheta(0)(3, 3) = exp(d4*2.0)*4.0;
    // ddtheta(0)(7, 7) = 2.0;
    // ddtheta(0)(8, 8) = 2.0;
    // ddtheta(0)(9, 9) = 2.0;
    // ddtheta(1)(0, 0) = t1*t5;
    // ddtheta(1)(0, 7) = t5;
    // ddtheta(1)(7, 0) = t5;
    // ddtheta(2)(1, 1) = t2*t6;
    // ddtheta(2)(1, 8) = t6;
    // ddtheta(2)(4, 7) = 1.0;
    // ddtheta(2)(7, 4) = 1.0;
    // ddtheta(2)(8, 1) = t6;
    // ddtheta(3)(2, 2) = t3*t7;
    // ddtheta(3)(2, 9) = t7;
    // ddtheta(3)(5, 8) = 1.0;
    // ddtheta(3)(6, 7) = 1.0;
    // ddtheta(3)(7, 6) = 1.0;
    // ddtheta(3)(8, 5) = 1.0;
    // ddtheta(3)(9, 2) = t7;
    // ddtheta(4)(1, 1) = t17;
    // ddtheta(4)(2, 2) = t18;
    // ddtheta(4)(4, 4) = 2.0;
    // ddtheta(4)(5, 5) = 2.0;
    // ddtheta(4)(6, 6) = 2.0;
    // ddtheta(5)(0, 0) = s12*t14;
    // ddtheta(5)(0, 4) = t14;
    // ddtheta(5)(4, 0) = t14;
    // ddtheta(6)(0, 0) = t16;
    // ddtheta(6)(2, 2) = t18;
    // ddtheta(6)(5, 5) = 2.0;
    // ddtheta(6)(6, 6) = 2.0;
    // ddtheta(7)(0, 0) = s13*t14;
    // ddtheta(7)(0, 6) = t14;
    // ddtheta(7)(6, 0) = t14;
    // ddtheta(8)(1, 1) = s23*t15;
    // ddtheta(8)(1, 5) = t15;
    // ddtheta(8)(4, 6) = -1.0;
    // ddtheta(8)(5, 1) = t15;
    // ddtheta(8)(6, 4) = -1.0;
    // ddtheta(9)(0, 0) = t16;
    // ddtheta(9)(1, 1) = t17;
    // ddtheta(9)(4, 4) = 2.0;

    Vec10d theta;

    const double alpha = z[0];
    const double d1 = z[1];
    const double d2 = z[2];
    const double d3 = z[3];
    const double s12 = z[4];
    const double s23 = z[5];
    const double s13 = z[6];
    const double t1 = z[7];
    const double t2 = z[8];
    const double t3 = z[9];

    const double exp_2alpha = std::exp(2 * alpha);
    const double exp_2d1 = std::exp(2 * d1);
    const double exp_2d2 = std::exp(2 * d2);
    const double exp_2d3 = std::exp(2 * d3);
    const double exp_d1 = std::exp(d1);
    const double exp_d2 = std::exp(d2);
    const double exp_d3 = std::exp(d3);

    theta[0] = 1;
    theta[1] = t1;
    theta[2] = t2;
    theta[3] = t3;
    theta[4] = std::pow(s23, 2) + std::pow(t2, 2) + std::pow(t3, 2) + std::pow(exp_d2, 2) + std::pow(exp_d3, 2);
    theta[5] = -s12 * exp_d2 - s13 * s23 - t1 * t2;
    theta[6] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(t1, 2) + std::pow(t3, 2) + std::pow(exp_d1, 2) + std::pow(exp_d3, 2);
    theta[7] = -s13 * exp_d3 - t1 * t3;
    theta[8] = -s23 * exp_d3 - t2 * t3;
    theta[9] = std::pow(s12, 2) + std::pow(s13, 2) + std::pow(s23, 2) + std::pow(t1, 2) + std::pow(t2, 2) + std::pow(exp_d1, 2) + std::pow(exp_d2, 2);
    const double exp_2_alpha = std::exp(2 * alpha);
    theta *= exp_2_alpha;

    dtheta.setZero();

    dtheta(0, 0) = 2 * exp_2alpha;

    dtheta(1, 0) = 2 * t1 * exp_2alpha;
    dtheta(1, 7) = exp_2alpha;

    dtheta(2, 0) = 2 * t2 * exp_2alpha;
    dtheta(2, 8) = exp_2alpha;

    dtheta(3, 0) = 2 * t3 * exp_2alpha;
    dtheta(3, 9) = exp_2alpha;

    dtheta(4, 0) = 2 * (std::pow(s23, 2) + std::pow(t2, 2) + std::pow(t3, 2) + exp_2d2 + exp_2d3) * exp_2alpha;
    dtheta(4, 2) = 2 * exp_2alpha * exp_2d2;
    dtheta(4, 3) = 2 * exp_2alpha * exp_2d3;
    dtheta(4, 5) = 2 * s23 * exp_2alpha;
    dtheta(4, 8) = 2 * t2 * exp_2alpha;
    dtheta(4, 9) = 2 * t3 * exp_2alpha;

    dtheta(5, 0) = -2 * (s12 * exp_d2 + s13 * s23 + t1 * t2) * exp_2alpha;
    dtheta(5, 2) = -s12 * exp_2alpha * exp_d2;
    dtheta(5, 4) = -exp_2alpha * exp_d2;
    dtheta(5, 5) = -s13 * exp_2alpha;
    dtheta(5, 6) = -s23 * exp_2alpha;
    dtheta(5, 7) = -t2 * exp_2alpha;
    dtheta(5, 8) = -t1 * exp_2alpha;

    dtheta(6, 0) = 2 * (std::pow(s12, 2) + std::pow(s13, 2) + std::pow(t1, 2) + std::pow(t3, 2) + exp_2d1 + exp_2d3) * exp_2alpha;
    dtheta(6, 1) = 2 * exp_2alpha * exp_2d1;
    dtheta(6, 3) = 2 * exp_2alpha * exp_2d3;
    dtheta(6, 4) = 2 * s12 * exp_2alpha;
    dtheta(6, 6) = 2 * s13 * exp_2alpha;
    dtheta(6, 7) = 2 * t1 * exp_2alpha;
    dtheta(6, 9) = 2 * t3 * exp_2alpha;

    dtheta(7, 0) = -2 * (s13 * exp_d3 + t1 * t3) * exp_2alpha;
    dtheta(7, 3) = -s13 * exp_2alpha * exp_d3;
    dtheta(7, 6) = -exp_2alpha * exp_d3;
    dtheta(7, 7) = -t3 * exp_2alpha;
    dtheta(7, 9) = -t1 * exp_2alpha;

    dtheta(8, 0) = -2 * (s23 * exp_d3 + t2 * t3) * exp_2alpha;
    dtheta(8, 3) = -s23 * exp_2alpha * exp_d3;
    dtheta(8, 5) = -exp_2alpha * exp_d3;
    dtheta(8, 8) = -t3 * exp_2alpha;
    dtheta(8, 9) = -t2 * exp_2alpha;

    dtheta(9, 0) = 2 * (std::pow(s12, 2) + std::pow(s13, 2) + std::pow(s23, 2) + std::pow(t1, 2) + std::pow(t2, 2) + exp_2d1 + exp_2d2) * exp_2alpha;
    dtheta(9, 1) = 2 * exp_2alpha * exp_2d1;
    dtheta(9, 2) = 2 * exp_2alpha * exp_2d2;
    dtheta(9, 4) = 2 * s12 * exp_2alpha;
    dtheta(9, 5) = 2 * s23 * exp_2alpha;
    dtheta(9, 6) = 2 * s13 * exp_2alpha;
    dtheta(9, 7) = 2 * t1 * exp_2alpha;
    dtheta(9, 8) = 2 * t2 * exp_2alpha;

    for (Index i = 0; i < 10; i++) {
        ddtheta(i).setZero();
    }

    double t5 = std::exp(d2);
    double t6 = std::exp(d3);
    double t7 = alpha*2.0;
    double t8 = d1*2.0;
    double t9 = d2*2.0;
    double t10 = d3*2.0;
    double t11 = s12*s12;
    double t12 = s13*s13;
    double t13 = s23*s23;
    double t14 = t1*t1;
    double t15 = t2*t2;
    double t16 = t3*t3;
    double t17 = std::exp(t7);
    double t18 = std::exp(t8);
    double t19 = std::exp(t9);
    double t20 = std::exp(t10);
    double t21 = t17*2.0;
    double t22 = t5*t17;
    double t23 = t6*t17;
    double t24 = -t17;
    double t26 = s12*t17*4.0;
    double t27 = s13*t17*4.0;
    double t29 = s23*t17*4.0;
    double t32 = t1*t17*4.0;
    double t34 = t2*t17*4.0;
    double t35 = t3*t17*4.0;
    double t38 = s13*t17*-2.0;
    double t39 = s23*t17*-2.0;
    double t40 = t1*t17*-2.0;
    double t41 = t2*t17*-2.0;
    double t42 = t3*t17*-2.0;
    double t50 = t17*t18*4.0;
    double t51 = t17*t19*4.0;
    double t52 = t17*t20*4.0;
    double t43 = -t22;
    double t44 = t22*-2.0;
    double t45 = -t23;
    double t46 = t23*-2.0;
    double t53 = s12*t44;
    double t54 = s13*t46;
    double t55 = s23*t46;
    ddtheta(0)(0, 0) = t17*4.0;
    ddtheta(1)(0, 0) = t32;
    ddtheta(1)(0, 7) = t21;
    ddtheta(1)(7, 0) = t21;
    ddtheta(2)(0, 0) = t34;
    ddtheta(2)(0, 8) = t21;
    ddtheta(2)(8, 0) = t21;
    ddtheta(3)(0, 0) = t35;
    ddtheta(3)(0, 9) = t21;
    ddtheta(3)(9, 0) = t21;
    ddtheta(4)(0, 0) = t17*(t13+t15+t16+t19+t20)*4.0;
    ddtheta(4)(0, 2) = t51;
    ddtheta(4)(0, 3) = t52;
    ddtheta(4)(0, 5) = t29;
    ddtheta(4)(0, 8) = t34;
    ddtheta(4)(0, 9) = t35;
    ddtheta(4)(2, 0) = t51;
    ddtheta(4)(2, 2) = t51;
    ddtheta(4)(3, 0) = t52;
    ddtheta(4)(3, 3) = t52;
    ddtheta(4)(5, 0) = t29;
    ddtheta(4)(5, 5) = t21;
    ddtheta(4)(8, 0) = t34;
    ddtheta(4)(8, 8) = t21;
    ddtheta(4)(9, 0) = t35;
    ddtheta(4)(9, 9) = t21;
    ddtheta(5)(0, 0) = t17*(s13*s23+s12*t5+t1*t2)*-4.0;
    ddtheta(5)(0, 2) = t53;
    ddtheta(5)(0, 4) = t44;
    ddtheta(5)(0, 5) = t38;
    ddtheta(5)(0, 6) = t39;
    ddtheta(5)(0, 7) = t41;
    ddtheta(5)(0, 8) = t40;
    ddtheta(5)(2, 0) = t53;
    ddtheta(5)(2, 2) = s12*t43;
    ddtheta(5)(2, 4) = t43;
    ddtheta(5)(4, 0) = t44;
    ddtheta(5)(4, 2) = t43;
    ddtheta(5)(5, 0) = t38;
    ddtheta(5)(5, 6) = t24;
    ddtheta(5)(6, 0) = t39;
    ddtheta(5)(6, 5) = t24;
    ddtheta(5)(7, 0) = t41;
    ddtheta(5)(7, 8) = t24;
    ddtheta(5)(8, 0) = t40;
    ddtheta(5)(8, 7) = t24;
    ddtheta(6)(0, 0) = t17*(t11+t12+t14+t16+t18+t20)*4.0;
    ddtheta(6)(0, 1) = t50;
    ddtheta(6)(0, 3) = t52;
    ddtheta(6)(0, 4) = t26;
    ddtheta(6)(0, 6) = t27;
    ddtheta(6)(0, 7) = t32;
    ddtheta(6)(0, 9) = t35;
    ddtheta(6)(1, 0) = t50;
    ddtheta(6)(1, 1) = t50;
    ddtheta(6)(3, 0) = t52;
    ddtheta(6)(3, 3) = t52;
    ddtheta(6)(4, 0) = t26;
    ddtheta(6)(4, 4) = t21;
    ddtheta(6)(6, 0) = t27;
    ddtheta(6)(6, 6) = t21;
    ddtheta(6)(7, 0) = t32;
    ddtheta(6)(7, 7) = t21;
    ddtheta(6)(9, 0) = t35;
    ddtheta(6)(9, 9) = t21;
    ddtheta(7)(0, 0) = t17*(s13*t6+t1*t3)*-4.0;
    ddtheta(7)(0, 3) = t54;
    ddtheta(7)(0, 6) = t46;
    ddtheta(7)(0, 7) = t42;
    ddtheta(7)(0, 9) = t40;
    ddtheta(7)(3, 0) = t54;
    ddtheta(7)(3, 3) = s13*t45;
    ddtheta(7)(3, 6) = t45;
    ddtheta(7)(6, 0) = t46;
    ddtheta(7)(6, 3) = t45;
    ddtheta(7)(7, 0) = t42;
    ddtheta(7)(7, 9) = t24;
    ddtheta(7)(9, 0) = t40;
    ddtheta(7)(9, 7) = t24;
    ddtheta(8)(0, 0) = t17*(s23*t6+t2*t3)*-4.0;
    ddtheta(8)(0, 3) = t55;
    ddtheta(8)(0, 5) = t46;
    ddtheta(8)(0, 8) = t42;
    ddtheta(8)(0, 9) = t41;
    ddtheta(8)(3, 0) = t55;
    ddtheta(8)(3, 3) = s23*t45;
    ddtheta(8)(3, 5) = t45;
    ddtheta(8)(5, 0) = t46;
    ddtheta(8)(5, 3) = t45;
    ddtheta(8)(8, 0) = t42;
    ddtheta(8)(8, 9) = t24;
    ddtheta(8)(9, 0) = t41;
    ddtheta(8)(9, 8) = t24;
    ddtheta(9)(0, 0) = t17*(t11+t12+t13+t14+t15+t18+t19)*4.0;
    ddtheta(9)(0, 1) = t50;
    ddtheta(9)(0, 2) = t51;
    ddtheta(9)(0, 4) = t26;
    ddtheta(9)(0, 5) = t29;
    ddtheta(9)(0, 6) = t27;
    ddtheta(9)(0, 7) = t32;
    ddtheta(9)(0, 8) = t34;
    ddtheta(9)(1, 0) = t50;
    ddtheta(9)(1, 1) = t50;
    ddtheta(9)(2, 0) = t51;
    ddtheta(9)(2, 2) = t51;
    ddtheta(9)(4, 0) = t26;
    ddtheta(9)(4, 4) = t21;
    ddtheta(9)(5, 0) = t29;
    ddtheta(9)(5, 5) = t21;
    ddtheta(9)(6, 0) = t27;
    ddtheta(9)(6, 6) = t21;
    ddtheta(9)(7, 0) = t32;
    ddtheta(9)(7, 7) = t21;
    ddtheta(9)(8, 0) = t34;
    ddtheta(9)(8, 8) = t21;

    return theta;
}

}; // namespace RAPTOR
