#include "EndEffectorMoment.h"

namespace RAPTOR {

// // constructor
// EndEffectorMoment::EndEffectorMoment()
// {
// }

// // destructors
// EndEffectorMoment::~EndEffectorMoment()
// {
// }

bool EndEffectorMoment::set_parameters(
    const std::shared_ptr<MatX>& A_full_input,
    const std::shared_ptr<VecX>& b_input,
    const std::shared_ptr<VecX>& inertia_parametersPtr_input,
    const bool include_offset_input
)
{ 
    enable_hessian = false;

    A = *A_full_input;
    b = *b_input;
    inertia_parametersPtr_ = inertia_parametersPtr_input;
    include_offset = include_offset_input;

    Nact = 7;

    // directly set up initial condition here
    int n = 10;  // End-effector parameters

    // initial guess for independent parameters
    // is just what is included in the optimised in the full parameters optimise
    x0 = VecX::Ones(n);
    std::cout << "end_effector parameters"<< x0.transpose() <<std::endl;

    return true;
}

bool EndEffectorMoment::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n =10;       // End-effector parameters 
    numVars= n;

    // number of inequality constraint
    numCons = 0;
    m = numCons;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool EndEffectorMoment::eval_f(
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
    // log-Cholesky parameterization
    double d1 = z[0];
    double d2 = z[1];
    double d3 = z[2];
    double d4 = z[3];
    double s12 = z[4];
    double s23 = z[5];
    double s13 = z[6];
    double t1 = z[7];
    double t2 = z[8];
    double t3 = z[9];

    Mat4 U;
    U << std::exp(d1), s12,      s13,      t1,
         0.0,     std::exp(d2),  s23,      t2,
         0.0,     0.0,      std::exp(d3),  t3,
         0.0,     0.0,      0.0,      std::exp(d4);

    // Compute LMI = U' * U
    Mat4 LMI = U.transpose() * U;

    // End-effector parameters
    VecX theta = VecX::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    // Update the inertia parameters
    inertia_parametersPtr_->tail(10) = theta;

    // Compute the ojective function
    VecX diff = A * (*inertia_parametersPtr_) - b;

    obj_value = 0.5 / b.size() * diff.transpose() * diff ;
    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool EndEffectorMoment::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);

    // log-Cholesky parameterization
    double d1 = z[0];
    double d2 = z[1];
    double d3 = z[2];
    double d4 = z[3];
    double s12 = z[4];
    double s23 = z[5];
    double s13 = z[6];
    double t1 = z[7];
    double t2 = z[8];
    double t3 = z[9];

    Mat4 U;
    U << std::exp(d1), s12,      s13,      t1,
         0.0,     std::exp(d2),  s23,      t2,
         0.0,     0.0,      std::exp(d3),  t3,
         0.0,     0.0,      0.0,      std::exp(d4);

    // Compute LMI = U' * U
    Mat4 LMI = U.transpose() * U;

    // End-effector parameters
    VecX theta = VecX::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    inertia_parametersPtr_->tail(10)= theta;
    VecX diff  = A * (*inertia_parametersPtr_) - b;

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

    MatX dtheta_dx = MatX::Zero(10, 10);

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
    grad_f_vec = (diff.transpose() * (A.rightCols(10) * dtheta_dx)) / b.size();

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i);
    }

    return true;
}
}; // namespace RAPTOR
