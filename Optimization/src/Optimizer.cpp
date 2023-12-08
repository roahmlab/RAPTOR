#include "Optimizer.h"

namespace IDTO {

// [TNLP_get_bounds_info]
// returns the variable bounds
bool Optimizer::get_bounds_info(
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
    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in get_bounds_info!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in get_bounds_info!");
    }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1e8;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1e8;
    }

    // compute bounds for all constraints
    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        constraintsPtrVec_[c]->compute_bounds();

        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g_l[iter] = constraintsPtrVec_[c]->g_lb(i);
            g_u[iter] = constraintsPtrVec_[c]->g_ub(i);
            iter++;
        }
    }

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool Optimizer::get_starting_point(
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

    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in get_starting_point!");
    }

    if (x0.size() != numVars) {
        throw std::runtime_error("You haven't specified an initial guess (x0) yet!");
    }

    for ( Index i = 0; i < n; i++ ) {
        x[i] = x0(i);
    }

    return true;
}
// [TNLP_get_starting_point]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool Optimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }

    // fill in a Eigen Vector instance of decision variables
    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    // auto start = std::chrono::high_resolution_clock::now();

    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // compute constraints
        try {
            constraintsPtrVec_[c]->compute(z, false);
        }
        catch (std::exception& e) {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("*** Error in eval_g!");
        }

        // fill in constraints
        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g[iter] = constraintsPtrVec_[c]->g(i);
            iter++;
        }
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << "g time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool Optimizer::eval_jac_g(
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
    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in eval_jac_g!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in eval_jac_g!");
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

        // auto start = std::chrono::high_resolution_clock::now();

        Index iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            // compute constraints
            try {
                constraintsPtrVec_[c]->compute(z, true);
            }
            catch (std::exception& e) {
                std::cout << e.what() << std::endl;
                throw std::runtime_error("*** Error in eval_jac_g!");
            }

            // fill in constraints
            for ( Index i = 0; i < constraintsPtrVec_[c]->pg_pz.rows(); i++ ) {
                for ( Index j = 0; j < constraintsPtrVec_[c]->pg_pz.cols(); j++ ) {
                    values[iter] = constraintsPtrVec_[c]->pg_pz(i, j);
                    iter++;
                }
            }
        }

        // auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "jac_g time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";
    }

    return true;
}
// [TNLP_eval_jac_g]


// [TNLP_eval_h]
//return the structure or values of the Hessian
bool Optimizer::eval_h(
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
void Optimizer::finalize_solution(
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
    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in finalize_solution!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in finalize_solution!");
    }

    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.

    // store the solution
    solution.resize(n);
    for( Index i = 0; i < n; i++ ) {
        solution(i) = x[i];
    }

    summarize_constraints(m, g);
}
// [TNLP_finalize_solution]

// [TNLP_summarize_constraints]
void Optimizer::summarize_constraints(
    Index                      m,
    const Number*              g
) 
{
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in summarize_constraints!");
    }

    std::cout << "Constraint violation report:" << std::endl;

    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // find where the maximum constraint violation is
        Number max_constr_violation = 0;
        Index max_constr_violation_id1 = 0;
        Index max_constr_violation_id2 = 0;
        for (Index i = 0; i < constraintsPtrVec_[c]->m; i++) {
            Number constr_violation = std::max(constraintsPtrVec_[c]->g_lb[i] - g[iter], g[iter] - constraintsPtrVec_[c]->g_ub[i]);
           
            if (constr_violation > max_constr_violation) {
                max_constr_violation_id1 = i;
                max_constr_violation_id2 = iter;
                max_constr_violation = constr_violation;
            }

            iter++;
        }

        // report constraint violation
        if (max_constr_violation > 0) {
            std::cout << constraintsNameVec_[c] << ": " << max_constr_violation << std::endl;
            std::cout << "    range: [" << constraintsPtrVec_[c]->g_lb[max_constr_violation_id1] 
                                        << ", " 
                                        << constraintsPtrVec_[c]->g_ub[max_constr_violation_id1] 
                      << "], value: "   << g[max_constr_violation_id2] << std::endl;
        }
    }
}
// [TNLP_summarize_constraints]

}; // namespace IDTO