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
    const int N_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    char stanceLeg, 
    const Transform& stance_foot_T_des_input
 ) 
{
    x0 = x0_input;

    fcPtr_ = std::make_unique<FourierCurves>(0.4, N_input, NUM_INDEPENDENT_JOINTS, Chebyshev, 6);

    dcidPtr_ = std::make_unique<DigitConstrainedInverseDynamics>(model_input, N_input, NUM_DEPENDENT_JOINTS, jtype_input, stanceLeg, stance_foot_T_des_input);

    assert(x0.size() == fcPtr_->varLength);

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
    n = fcPtr_->varLength;

    // number of inequality constraint
    m = fcPtr_->Nact * 3;

    VecX z0(n);
    for ( Index i = 0; i < n; i++ ) {
        z0(i) = x0[i];
    }

    fcPtr_->compute(z0, true);

    nnz_jac_g = fcPtr_->pq_pz(20).nonZeros() + 
                fcPtr_->pq_d_pz(20).nonZeros() + 
                fcPtr_->pq_dd_pz(20).nonZeros();

    // printf("nnz_jac_g: %d\n", nnz_jac_g);
    // printf("nnz_jac_g 1: %ld\n", fcPtr_->pq_pz(20).nonZeros());
    // printf("nnz_jac_g 2: %ld\n", fcPtr_->pq_d_pz(20).nonZeros());
    // printf("nnz_jac_g 3: %ld\n", fcPtr_->pq_dd_pz(20).nonZeros());

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
    if(n != fcPtr_->varLength){
        throw std::runtime_error("*** Error wrong value of n in get_bounds_info!");
    }
    // if(m != fcPtr_->Nact * 3){
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

    for( Index i = 0; i < m; i++ ) {
        g_l[i] = -1e8;
        g_u[i] = 1e8;
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
    if(init_x == false || init_z == true || init_lambda == true){
        throw std::runtime_error("*** Error wrong value of init in get_starting_point!");
    }

    if(n != fcPtr_->varLength){
        throw std::runtime_error("*** Error wrong value of n in get_starting_point!");
    }

    for( Index i = 0; i < n; i++ ) {
        // initialize to zero
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
    if(n != fcPtr_->varLength){
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
    if(n != fcPtr_->varLength){
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
    if(n != fcPtr_->varLength){
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if(m != fcPtr_->Nact * 3){
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }

    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];

        // cout << x[i] << endl;
    }
    // cout << endl;

    // if (new_x) {
        fcPtr_->compute(z, false);
    // }

    for ( Index i = 0; i < fcPtr_->Nact; i++ ) {
        g[i] = fcPtr_->q(20)(i);
        g[i + fcPtr_->Nact] = fcPtr_->q_d(20)(i);
        g[i + 2 * fcPtr_->Nact] = fcPtr_->q_dd(20)(i);
    }

    // cout << fcPtr_->q(20).transpose() << endl;
    // cout << fcPtr_->q_d(20).transpose() << endl;
    // cout << fcPtr_->q_dd(20).transpose() << endl;
    // cout << endl;

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
    if(n != fcPtr_->varLength){
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if(m != fcPtr_->Nact * 3){
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }
        
    if( values == NULL ) {
        fcPtr_->compute(x0, true);

        // return the structure of the Jacobian
        Index iter = 0;

        const SpaMatX& pq_pz = fcPtr_->pq_pz(20);
        for (Index i = 0; i < pq_pz.outerSize(); i++) {
            Index k_start = pq_pz.outerIndexPtr()[i];
            Index k_end   = pq_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                Index j = pq_pz.innerIndexPtr()[k];
                // double v = pq_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                iRow[iter] = j;
                jCol[iter] = i;
                iter++;
            }
        }

        const SpaMatX& pq_d_pz = fcPtr_->pq_d_pz(20);
        for (Index i = 0; i < pq_d_pz.outerSize(); i++) {
            Index k_start = pq_d_pz.outerIndexPtr()[i];
            Index k_end   = pq_d_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                Index j = pq_d_pz.innerIndexPtr()[k];
                // double v = pq_d_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                iRow[iter] = j + fcPtr_->Nact;
                jCol[iter] = i;
                iter++;
            }
        }

        const SpaMatX& pq_dd_pz = fcPtr_->pq_dd_pz(20);
        for (Index i = 0; i < pq_dd_pz.outerSize(); i++) {
            Index k_start = pq_dd_pz.outerIndexPtr()[i];
            Index k_end   = pq_dd_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                Index j = pq_dd_pz.innerIndexPtr()[k];
                // double v = pq_dd_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                iRow[iter] = j + 2 * fcPtr_->Nact;
                jCol[iter] = i;
                iter++;
            }
        }

        // cout << iter << endl;
        // for (Index i = 0; i < iter; i++) {
        //     cout << iRow[i] << " " << jCol[i] << endl;
        // }
        // cout << endl;
    }
    else {
        VecX z(n);
        for ( Index i = 0; i < n; i++ ) {
            z(i) = x[i];
        }

        fcPtr_->compute(z, true);

        Index iter = 0;

        const SpaMatX& pq_pz = fcPtr_->pq_pz(20);
        for (Index i = 0; i < pq_pz.outerSize(); i++) {
            Index k_start = pq_pz.outerIndexPtr()[i];
            Index k_end   = pq_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                // Index j = pq_pz.innerIndexPtr()[k];
                double v = pq_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                values[iter] = v;
                iter++;
            }
        }

        const SpaMatX& pq_d_pz = fcPtr_->pq_d_pz(20);
        for (Index i = 0; i < pq_d_pz.outerSize(); i++) {
            Index k_start = pq_d_pz.outerIndexPtr()[i];
            Index k_end   = pq_d_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                // Index j = pq_d_pz.innerIndexPtr()[k];
                double v = pq_d_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                values[iter] = v;
                iter++;
            }
        }

        const SpaMatX& pq_dd_pz = fcPtr_->pq_dd_pz(20);
        for (Index i = 0; i < pq_dd_pz.outerSize(); i++) {
            Index k_start = pq_dd_pz.outerIndexPtr()[i];
            Index k_end   = pq_dd_pz.outerIndexPtr()[i+1];

            for (Index k = k_start; k < k_end; k++) {
                // Index j = pq_dd_pz.innerIndexPtr()[k];
                double v = pq_dd_pz.valuePtr()[k];
                // v is value of the element at position (j,i)

                values[iter] = v;
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