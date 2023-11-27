#include "DigitSingleStepOptimizer.h"

namespace IDTO {
namespace Digit {

using std::cout;
using std::endl;

// // constructor
// DigitSingleStepOptimizer::DigitSingleStepOptimizer()
// {
// }


// // destructor
// DigitSingleStepOptimizer::~DigitSingleStepOptimizer()
// {
// }


bool DigitSingleStepOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    char stanceLeg, 
    const Transform& stance_foot_T_des_input
 ) 
{
    x0 = x0_input;

    fcPtr_ = std::make_shared<FourierCurves>(T_input, 
                                             N_input, 
                                             NUM_INDEPENDENT_JOINTS, 
                                             Chebyshev, 
                                             degree_input);

    dcidPtr_ = std::make_shared<DigitConstrainedInverseDynamics>(model_input, 
                                                                 N_input, 
                                                                 NUM_DEPENDENT_JOINTS, 
                                                                 jtype_input, 
                                                                 stanceLeg, 
                                                                 stance_foot_T_des_input);                                                          

    // convert to their base class pointers
    trajPtr_ = fcPtr_;
    dcPtr_ = dcidPtr_->dynamicsConstraintsPtr_;

    VecX JOINT_LIMITS_LOWER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_LOWER_VEC(i) = deg2rad(JOINT_LIMITS_LOWER[i]);
    }

    VecX JOINT_LIMITS_UPPER_VEC(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        JOINT_LIMITS_UPPER_VEC(i) = deg2rad(JOINT_LIMITS_UPPER[i]);
    }

    constraints_cjlPtr_ = std::make_unique<ConstrainedJointLimits>(trajPtr_, 
                                                                   dcPtr_, 
                                                                   JOINT_LIMITS_LOWER_VEC, 
                                                                   JOINT_LIMITS_UPPER_VEC);

    assert(x0.size() == trajPtr_->varLength);

    return true;
}
// [TNLP_set_parameters]

bool DigitSingleStepOptimizer::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // number of decision variables
    n = trajPtr_->varLength;

    // number of inequality constraint
    m = constraints_cjlPtr_->m;

    VecX z0(n);
    for ( Index i = 0; i < n; i++ ) {
        z0(i) = x0[i];
    }

    nnz_jac_g = n * m;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool DigitSingleStepOptimizer::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    if (n != trajPtr_->varLength) {
        throw std::runtime_error("*** Error wrong value of n in get_bounds_info!");
    }
    // if (m != constraints_cjlPtr_->m) {
    //     throw std::runtime_error("*** Error wrong value of m in get_bounds_info!");
    // }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1e8;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1e8;
    }

    // compute bounds for all constraints
    constraints_cjlPtr_->compute_bounds();

    // fill in bounds for all constraints
    Index iter = 0;
    for( Index i = 0; i < constraints_cjlPtr_->m; i++ ) {
        g_l[iter] = constraints_cjlPtr_->g_lb(i);
        g_u[iter] = constraints_cjlPtr_->g_ub(i);
        iter++;
    }

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool DigitSingleStepOptimizer::get_starting_point(
    Index   n,
    bool    init_x,
    Number* x,
    bool    init_z,
    Number* z_L,
    Number* z_U,
    Index   m,
    bool    init_lambda,
    Number* lambda
)
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    if (init_x == false || init_z == true || init_lambda == true) {
        throw std::runtime_error("*** Error wrong value of init in get_starting_point!");
    }

    if (n != trajPtr_->varLength) {
        throw std::runtime_error("*** Error wrong value of n in get_starting_point!");
    }

    for ( Index i = 0; i < n; i++ ) {
        x[i] = x0(i);
    }

    return true;
}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool DigitSingleStepOptimizer::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != trajPtr_->varLength){
       throw std::runtime_error("*** Error wrong value of n in eval_f!");
    }

    obj_value = 0;

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DigitSingleStepOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != trajPtr_->varLength){
       throw std::runtime_error("*** Error wrong value of n in eval_f!");
    }

    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] = 0;
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool DigitSingleStepOptimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if (n != trajPtr_->varLength) {
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if (m != constraints_cjlPtr_->m) {
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }

    // fill in a Eigen Vector instance of decision variables
    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    // compute all constraints
    try {
        constraints_cjlPtr_->compute(z, false);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        throw std::runtime_error("*** Error in eval_g!");
    }

    // fill in all constraints
    Index iter = 0;
    for ( Index i = 0; i < constraints_cjlPtr_->m; i++ ) {
        g[iter] = constraints_cjlPtr_->g(i);
        iter++;
    }

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool DigitSingleStepOptimizer::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    if (n != trajPtr_->varLength) {
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if (m != constraints_cjlPtr_->m) {
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }
        
    if( values == NULL ) {
        // return the structure of the Jacobian
        // this particular Jacobian is dense
        for(Index i = 0; i < m; i++){
            for(Index j = 0; j < n; j++){
                iRow[i * n + j] = i;
                jCol[i * n + j] = j;
            }
        }
    }
    else {
        // fill in a Eigen Vector instance of decision variables
        VecX z(n);
        for ( Index i = 0; i < n; i++ ) {
            z(i) = x[i];
        }

        // compute all constraints
        try {
            constraints_cjlPtr_->compute(z, true);
        }
        catch (std::exception& e) {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("*** Error in eval_jac_g!");
        }

        // quick sanity check
        if (constraints_cjlPtr_->pg_pz.rows() != constraints_cjlPtr_->m) {
            throw std::runtime_error("Error wrong value of pg_pz.rows() in eval_jac_g!");
        }
        if (constraints_cjlPtr_->pg_pz.cols() != n) {
            throw std::runtime_error("Error wrong value of pg_pz.cols() in eval_jac_g!");
        }
        
        // fill in all constraints
        Index iter = 0;

        for ( Index i = 0; i < constraints_cjlPtr_->pg_pz.rows(); i++ ) {
            for ( Index j = 0; j < constraints_cjlPtr_->pg_pz.cols(); j++ ) {
                values[iter] = constraints_cjlPtr_->pg_pz(i, j);
                iter++;
            }
        }
    }

    return true;
}
// [TNLP_eval_jac_g]


// [TNLP_eval_h]
//return the structure or values of the Hessian
bool DigitSingleStepOptimizer::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    return false;
}
// [TNLP_eval_h]


// [TNLP_finalize_solution]
void DigitSingleStepOptimizer::finalize_solution(
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
)
{
    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.

    // store the solution
    // for( Index i = 0; i < n; i++ ) {
    //     solution[i] = (double)x[i];
    // }

    
}
// [TNLP_finalize_solution]

}; // namespace Digit
}; // namespace IDTO