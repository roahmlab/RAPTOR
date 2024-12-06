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
    const VecXd offset_input
)
{ 
    enable_hessian = false;

    // parse the robot model
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(*modelPtr_);

    phi = VecXd::Zero(10 * modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
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
    trajPtr_ = std::make_shared<TrajectoryData>(filename_input);

    // this trajectory is to compute gravity regressors,
    // so set velocity to 0 while acceleration is already 0 in TrajectoryData
    trajPtr2_ = std::make_shared<TrajectoryData>(filename_input);
    for (int i = 0; i < trajPtr2_->N; i++) {
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
        throw std::invalid_argument("0 segments");
    }

    // trajectoryData does not require any decision variable, 
    // so simply put an empty Eigen vector as placeholder here
    mrPtr_->compute(VecXd::Zero(1), false);
    ridPtr_->compute(VecXd::Zero(1), false);

    // now compute regression elements A and b
    // which are essentially combination of different dynamic regressors
    A.resize(modelPtr_->nv * num_segment, 10 * modelPtr_->nv);
    b.resize(modelPtr_->nv * num_segment);

    int i = 0;
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
    sensor_noise = sensor_noise_input;

    if ((IntervalHelper::getCenter(sensor_noise.position_error) == 0 &&
         IntervalHelper::getRadius(sensor_noise.position_error) == 0) &&
        (IntervalHelper::getCenter(sensor_noise.velocity_error) == 0 &&
         IntervalHelper::getRadius(sensor_noise.velocity_error) == 0) &&
        (IntervalHelper::getCenter(sensor_noise.acceleration_error) == 0 &&
         IntervalHelper::getRadius(sensor_noise.acceleration_error) == 0)) {
        std::cout << "No sensor noise is added" << std::endl;
    }
    else {
        std::cout << "Sensor noise is added" << std::endl;
        std::cout << "Position error: " << sensor_noise.position_error.lower() 
                                 << ' ' << sensor_noise.position_error.upper() << std::endl;
        std::cout << "Velocity error: " << sensor_noise.velocity_error.lower() 
                                 << ' ' << sensor_noise.velocity_error.upper() << std::endl;

        mrIntPtr_ = std::make_shared<IntervalMomentumRegressor>(*modelPtr_, trajPtr_, sensor_noise);
        ridIntPtr_ = std::make_shared<IntervalRegressorInverseDynamics>(*modelPtr_, trajPtr2_, sensor_noise);
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

    // log-Cholesky parameterization
    const double d1 = z[0];
    const double d2 = z[1];
    const double d3 = z[2];
    const double d4 = z[3];
    const double s12 = z[4];
    const double s23 = z[5];
    const double s13 = z[6];
    const double t1 = z[7];
    const double t2 = z[8];
    const double t3 = z[9];

    Mat4d U;
    U << std::exp(d1), s12,          s13,          t1,
         0.0,          std::exp(d2), s23,          t2,
         0.0,          0.0,          std::exp(d3), t3,
         0.0,          0.0,          0.0,          std::exp(d4);

    // Compute LMI = U' * U
    Mat4d LMI = U.transpose() * U;

    // End-effector parameters
    VecXd theta = VecXd::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    // Update the inertia parameters
    phi.tail(10) = theta;

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
    VecXd z = Utils::initializeEigenVectorFromArray(x, n);
    VecXd grad_f_vec = VecXd::Zero(n);

    // log-Cholesky parameterization
    const double d1 = z[0];
    const double d2 = z[1];
    const double d3 = z[2];
    const double d4 = z[3];
    const double s12 = z[4];
    const double s23 = z[5];
    const double s13 = z[6];
    const double t1 = z[7];
    const double t2 = z[8];
    const double t3 = z[9];

    Mat4d U;
    U << std::exp(d1), s12,      s13,      t1,
         0.0,     std::exp(d2),  s23,      t2,
         0.0,     0.0,      std::exp(d3),  t3,
         0.0,     0.0,      0.0,      std::exp(d4);

    // Compute LMI = U' * U
    Mat4d LMI = U.transpose() * U;

    // End-effector parameters
    VecXd theta = VecXd::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    phi.tail(10)= theta;
    VecXd diff = A * phi - b;

    // Compute the gradient of dtheta_ dx by matlab symbolic toolbox
    double t5 = std::exp(d1);
    double t6 = std::exp(d2);
    double t7 = std::exp(d3);
    double t8 = d1*2.0;
    double t9 = d2*2.0;
    double t10 = d3*2.0;
    double t11 = s12*2.0;
    double t12 = s13*2.0;
    double t13 = s23*2.0;
    double t14 = std::exp(t8);
    double t15 = std::exp(t9);
    double t16 = std::exp(t10);
    double t17 = -t5;
    double t18 = t14*2.0;
    double t19 = t15*2.0;
    double t20 = t16*2.0;

    MatXd dtheta_dx = MatXd::Zero(10, 10);
    dtheta_dx(0, 3) = std::exp(d4 * 2.0) * 2.0;
    dtheta_dx(0, 7) = t1 * 2.0;
    dtheta_dx(0, 8) = t2 * 2.0;
    dtheta_dx(0, 9) = t3 * 2.0;
    dtheta_dx(1, 0) = t1 * t5;
    dtheta_dx(1, 7) = t5;
    dtheta_dx(2, 1) = t2 * t6;
    dtheta_dx(2, 4) = t1;
    dtheta_dx(2, 7) = s12;
    dtheta_dx(2, 8) = t6;
    dtheta_dx(3, 2) = t3 * t7;
    dtheta_dx(3, 5) = t2;
    dtheta_dx(3, 6) = t1;
    dtheta_dx(3, 7) = s13;
    dtheta_dx(3, 8) = s23;
    dtheta_dx(3, 9) = t7;
    dtheta_dx(4, 1) = t19;
    dtheta_dx(4, 2) = t20;
    dtheta_dx(4, 4) = t11;
    dtheta_dx(4, 5) = t13;
    dtheta_dx(4, 6) = t12;
    dtheta_dx(5, 0) = s12 * t17;
    dtheta_dx(5, 4) = t17;
    dtheta_dx(6, 0) = t18;
    dtheta_dx(6, 2) = t20;
    dtheta_dx(6, 5) = t13;
    dtheta_dx(6, 6) = t12;
    dtheta_dx(7, 0) = s13 * t17;
    dtheta_dx(7, 6) = t17;
    dtheta_dx(8, 1) = -s23 * t6;
    dtheta_dx(8, 4) = -s13;
    dtheta_dx(8, 5) = -t6;
    dtheta_dx(8, 6) = -s12;
    dtheta_dx(9, 0) = t18;
    dtheta_dx(9, 1) = t19;
    dtheta_dx(9, 4) = t11;

    // Compute the gradient
    grad_f_vec = (diff.transpose() * A.rightCols(10) * dtheta_dx) / b.size();

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i);
    }

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

    if (mrIntPtr_ == nullptr || ridIntPtr_ == nullptr) {
        return;
    }

    std::cout << "Performing error analysis" << std::endl;

    // now compute interval version of regression elements A and b
    // which are essentially combination of different dynamic regressors
    mrIntPtr_->compute(VecXd::Zero(1), true);
    ridIntPtr_->compute(VecXd::Zero(1), true);

    Aint.resize(modelPtr_->nv * num_segment, 10 * modelPtr_->nv);
    bint.resize(modelPtr_->nv * num_segment);

    Aint2.resize(modelPtr_->nv * num_segment, 10 * modelPtr_->nv);

    int i = 0;
    #pragma omp parallel for shared(modelPtr_, trajPtr_,  mrIntPtr_, ridIntPtr_, Aint, bint) private(i) schedule(dynamic, 1)
    for (i = 0; i < 256 - H - 1; i += H) {
        int seg_start = i;
        int seg_end = seg_start + H;

        const MatXInt& Y_Hqd_1 = mrIntPtr_->Y.middleRows(seg_start * modelPtr_->nv, modelPtr_->nv);
        const MatXInt& Y_Hqd_2 = mrIntPtr_->Y.middleRows(seg_end * modelPtr_->nv, modelPtr_->nv);

        MatXInt int_Y_CTqd_g = MatXInt::Zero(modelPtr_->nv, 10 * modelPtr_->nv);
        VecXInt int_ctrl = VecXInt::Zero(modelPtr_->nv);

        for (int j = seg_start; j < seg_end; j++) {
            double dt = trajPtr_->tspan(j + 1) - trajPtr_->tspan(j);

            MatXInt Y_CTqd_i = mrIntPtr_->Y_CTv.middleRows(j * modelPtr_->nv, modelPtr_->nv);
            MatXInt Yg_i = ridIntPtr_->Y.middleRows(j * modelPtr_->nv, modelPtr_->nv);

            int_Y_CTqd_g += (Y_CTqd_i - Yg_i) * dt;

            // Note that here trajPtr_->q_dd stores the applied torque
            for (int k = 0; k < modelPtr_->nv; k++) {
                const Interval q_d_int = trajPtr_->q_d(j)(k) + sensor_noise.velocity_error;
                const Interval q_dd_int = trajPtr_->q_dd(j)(k) + sensor_noise.acceleration_error;
                
                int_ctrl(k) += (q_dd_int - 
                                modelPtr_->friction(k) * Utils::sign(trajPtr_->q_d(j)(k)) - 
                                modelPtr_->damping(k) * q_d_int - 
                                offset(k)) * dt;
            }
        }

        for (int k = 0; k < modelPtr_->nv; k++) {
            const Interval q_d_seg_start_int = trajPtr_->q_d(seg_start)(k) + sensor_noise.velocity_error;
            const Interval q_d_seg_end_int = trajPtr_->q_d(seg_end)(k) + sensor_noise.velocity_error;
            int_ctrl(k) -= modelPtr_->armature(k) * (q_d_seg_end_int - q_d_seg_start_int);
        }
        
        int s = i / H;
        Aint.middleRows(s * modelPtr_->nv, modelPtr_->nv) = (Y_Hqd_2 - Y_Hqd_1) - int_Y_CTqd_g;
        bint.segment(s * modelPtr_->nv, modelPtr_->nv) = int_ctrl;
    }

    // fast check
    int index = 0;
    for (int i = 60; i < Aint.rows(); i++) {
        for (int j = 0; j < Aint.cols(); j++) {
            // if (A(i, j) > Aint(i, j).upper() || 
            //     A(i, j) < Aint(i, j).lower()) {
            //     std::cerr << i << ' ' << j << std::endl;
            //     std::cerr << A(i, j) << " [" << Aint(i, j).lower() << ", " << Aint(i, j).upper() << "]" << std::endl;
            //     throw std::runtime_error("Wrong computation of Aint");
            // }

            // if (index <= 1000) {
            //     std::cout << " [" << Aint(i, j).lower() << ", " << Aint(i, j).upper() << "]" << std::endl;
            //     index++;
            // }

            double center = IntervalHelper::getCenter(Aint(i, j));
            double radius = IntervalHelper::getRadius(Aint(i, j));

            if (center != 0) {
                std::cout << center << ' ' << radius << ' ' << radius / center << std::endl;
            }
        }
    }

    // // for (int i = 0; i < bint.size(); i++) {
    // //     if (b(i) > bint(i).upper() || 
    // //         b(i) < bint(i).lower()) {
    // //         std::cerr << i << std::endl;
    // //         std::cerr << b(i) << " [" << bint(i).lower() << ", " << bint(i).upper() << "]" << std::endl;
    // //         throw std::runtime_error("Wrong computation of bint");
    // //     }
    // // }

    MatX AAT = A.rightCols(10).transpose()  * A.rightCols(10);
    MatXInt AATint = intervalMatrixMultiply(Aint.rightCols(10).transpose(), Aint.rightCols(10));
    for (int i = 0; i < AATint.rows(); i++) {
        for (int j = 0; j < AATint.cols(); j++) {
            std::cout << " [" << AATint(i, j).lower() << ", " << AATint(i, j).upper() << "]  ";
        }
        std::cout << std::endl;
    }

    MatX AAT_perturb(10, 10);
    for (int i = 0; i < AATint.rows(); i++) {
        for (int j = 0; j < AATint.cols(); j++) {
            AAT_perturb(i, j) = 0.5 * IntervalHelper::getRadius(AATint(i, j));
        }
    }

    std::cout << AAT_perturb << std::endl;
}

}; // namespace RAPTOR
