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
    const VecX offset_input
)
{ 
    enable_hessian = true;

    // parse the robot model
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(*modelPtr_);

    phi = VecX::Zero(10 * modelPtr_->nv);
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
        offset = VecX::Zero(modelPtr_->nv);
    }

    // simply give 0 as initial guess
    x0 = VecX::Zero(10);

    return true;
}

void EndEffectorParametersIdentification::add_trajectory_file(
    const std::shared_ptr<TrajectoryData>& trajectory_input,
    const MatX& acceleration_input) {
    trajPtrs_.push_back(trajectory_input);

    if (trajPtrs_.back()->Nact != modelPtr_->nv) {
        throw std::invalid_argument("The number of active joints in the trajectory does not match the robot model.");
    }

    accelerations_.push_back(acceleration_input);

    if (accelerations_.back().rows() != trajPtrs_.back()->N) {
        throw std::invalid_argument("The number of rows in the acceleration file does not match the trajectory.");
    }

    if (accelerations_.back().cols() != modelPtr_->nv) {
        throw std::invalid_argument("The number of columns in the acceleration file does not match the robot model.");
    }

    initialize_regressors(trajPtrs_.back(),
                          accelerations_.back());

    // combine all the regressors
    num_segments.clear();
    total_num_segments = 0;
    for (Index i = 0; i < Aseg.size(); i++) {
        num_segments.push_back(Aseg[i].rows() / modelPtr_->nv);
        total_num_segments += num_segments.back();
    }
}

void EndEffectorParametersIdentification::add_trajectory_file(
    const std::string hardware_trajectory_filename_input,
    const std::string acceleration_filename_input,
    const TimeFormat time_format,
    const int downsample_rate) {
    trajPtrs_.push_back(std::make_shared<TrajectoryData>(hardware_trajectory_filename_input, 
                                                         SensorNoiseInfo(modelPtr_->nv),
                                                         time_format,
                                                         downsample_rate));

    if (trajPtrs_.back()->Nact != modelPtr_->nv) {
        throw std::invalid_argument("The number of active joints in the trajectory does not match the robot model.");
    }

    accelerations_.push_back(Utils::initializeEigenMatrixFromFile(acceleration_filename_input));

    if (accelerations_.back().rows() != trajPtrs_.back()->N) {
        throw std::invalid_argument("The number of rows in the acceleration file does not match the trajectory.");
    }

    if (accelerations_.back().cols() != modelPtr_->nv) {
        throw std::invalid_argument("The number of columns in the acceleration file does not match the robot model.");
    }

    trajectoryFilenames_.push_back(hardware_trajectory_filename_input);

    initialize_regressors(trajPtrs_.back(),
                          accelerations_.back());

    // combine all the regressors
    num_segments.clear();
    total_num_segments = 0;
    for (Index i = 0; i < Aseg.size(); i++) {
        num_segments.push_back(Aseg[i].rows() / modelPtr_->nv);
        total_num_segments += num_segments.back();
    }
}

void EndEffectorParametersIdentification::initialize_regressors(const std::shared_ptr<TrajectoryData>& trajPtr_,
                                                                const MatX& acceleration) {
    MatX A_seg_i = MatX::Zero(modelPtr_->nv * trajPtr_->N, 10 * modelPtr_->nv);
    VecX b_seg_i = VecX::Zero(modelPtr_->nv * trajPtr_->N);

    Index i = 0;
    // #pragma omp parallel for shared(modelPtr_, dataPtr_, trajPtr_, A_seg_i, b_seg_i) private(i) schedule(dynamic, 1)
    for (i = 0; i < trajPtr_->N; i++) {
        // the following might be really confusing, but trajPtr_->q_dd stores the torque data from the hardware
        // and acceleration stores the actual acceleration data from a separate file
        const VecX& q = trajPtr_->q(i);
        const VecX& q_d = trajPtr_->q_d(i);
        const VecX& q_dd = acceleration.row(i);
        const VecX& tau = trajPtr_->q_dd(i);

        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, q, q_d, q_dd);

        A_seg_i.middleRows(i * modelPtr_->nv, modelPtr_->nv) = dataPtr_->jointTorqueRegressor;
        b_seg_i.segment(i * modelPtr_->nv, modelPtr_->nv) = tau -
                                                            modelPtr_->friction.cwiseProduct(q_d.cwiseSign()) -
                                                            modelPtr_->damping.cwiseProduct(q_d) -
                                                            modelPtr_->armature.cwiseProduct(q_dd) -
                                                            offset;


    }

    Aseg.push_back(A_seg_i);
    bseg.push_back(b_seg_i);
}

void EndEffectorParametersIdentification::reset() {
    Optimizer::reset();

    trajPtrs_.clear();
    accelerations_.clear();
    trajectoryFilenames_.clear();

    Aseg.clear();
    bseg.clear();

    num_segments.clear();

    std::cout << "End effector identification reset, all trajectory data removed." << std::endl;
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

    // number of constraints
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

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    phi.tail(10) = z_to_theta(z);

    // Compute the ojective function
    obj_value = 0;
    for (Index i = 0; i < Aseg.size(); i++) {
        const VecX diff = Aseg[i] * phi - bseg[i];
        obj_value += 0.5 * diff.dot(diff);
    }

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

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    Mat10 dtheta;
    phi.tail(10) = d_z_to_theta(z, dtheta);

    // Compute the gradient
    VecX grad_f_vec = VecX::Zero(n);
    for (Index i = 0; i < Aseg.size(); i++) {
        const VecX diff = Aseg[i] * phi - bseg[i];
        grad_f_vec += diff.transpose() * Aseg[i].rightCols(10) * dtheta;
    }

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

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    Mat10 dtheta;
    Eigen::Array<Mat10, 1, 10> ddtheta;
    phi.tail(10) = dd_z_to_theta(z, dtheta, ddtheta);

    // Compute the Hessian
    hess_f = MatX::Zero(n, n);
    for (Index i = 0; i < Aseg.size(); i++) {
        const VecX diff = Aseg[i] * phi - bseg[i];
        MatX temp1 = Aseg[i].rightCols(10) * dtheta;
        hess_f += temp1.transpose() * temp1;

        MatX temp2 = diff.transpose() * Aseg[i].rightCols(10);
        for (Index j = 0; j < n; j++) {
            hess_f += temp2(j) * ddtheta(j);
        }
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

    theta_solution = z_to_theta(solution);
}

Eigen::Vector<double, 10> EndEffectorParametersIdentification::z_to_theta(const VecX& z) {
    Vec10 theta;
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
    const VecX& z,
    Mat10& dtheta) {
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

    Vec10 theta;
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
    const VecX& z,
    Mat10& dtheta,
    Eigen::Array<Mat10, 1, 10>& ddtheta) {
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

    Vec10 theta;
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
