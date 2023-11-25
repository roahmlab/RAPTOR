#include "DigitSingleStepOptimizer.h"

namespace IDTO {
namespace Digit {

// // constructor
// DigitSingleStepOptimizer::DigitSingleStepOptimizer()
// {
// }


// // destructor
// DigitSingleStepOptimizer::~DigitSingleStepOptimizer()
// {
//     delete[] g_copy;
// }


bool DigitSingleStepOptimizer::set_parameters(
 ) 
{
    

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
    // // The problem described NUM_FACTORS variables, x[NUM_FACTORS] through x[NUM_FACTORS] for each joint
    // n = NUM_FACTORS;

    // // number of inequality constraint
    // m = constraint_number;

    // nnz_jac_g = m * n;

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
    // if(n != NUM_FACTORS){
    //     WARNING_PRINT("*** Error wrong value of n in get_bounds_info!");
    // }
    // if(m != constraint_number){
    //     WARNING_PRINT("*** Error wrong value of m in get_bounds_info!");
    // }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1.0;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1.0;
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

    // if(n != NUM_FACTORS){
    //     WARNING_PRINT("*** Error wrong value of n in get_starting_point!");
    // }

    for( Index i = 0; i < n; i++ ) {
        // initialize to zero
        x[i] = 0.0;
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
    // if(n != NUM_FACTORS){
    //    WARNING_PRINT("*** Error wrong value of n in eval_f!");
    // }

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
    // if(n != NUM_FACTORS){
    //     WARNING_PRINT("*** Error wrong value of n in eval_grad_f!");
    // }

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
    // if(n != NUM_FACTORS){
    //     WARNING_PRINT("*** Error wrong value of n in eval_g!");
    // }
    // if(m != constraint_number){
    //     WARNING_PRINT("*** Error wrong value of m in eval_g!");
    // }

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
    // if(n != NUM_FACTORS){
    //     WARNING_PRINT("*** Error wrong value of n in eval_g!");
    // }
    // if(m != constraint_number){
    //     WARNING_PRINT("*** Error wrong value of m in eval_g!");
    // }
        
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