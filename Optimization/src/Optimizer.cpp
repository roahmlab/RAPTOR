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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_bounds_info!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in get_bounds_info!");
    }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1e19;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1e19;
    }

    if (constraintsPtrVec_.size() != constraintsNameVec_.size()) {
        THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ and constraintsNameVec_ have different sizes!");
    }

    // compute bounds for all constraints
    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        try {
            constraintsPtrVec_[c]->compute_bounds();
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in get_bounds_info! Check previous error message.");
        }

        if (constraintsPtrVec_[c]->m != constraintsPtrVec_[c]->g_lb.size() || 
            constraintsPtrVec_[c]->m != constraintsPtrVec_[c]->g_ub.size()) {
            THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ have different sizes!");
        }

        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g_l[iter] = constraintsPtrVec_[c]->g_lb(i);
            g_u[iter] = constraintsPtrVec_[c]->g_ub(i);
            iter++;
        }
    }

    // report constraints distribution
    std::cout << "Dimension of each constraints and their locations: \n";
    iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        std::cout << constraintsNameVec_[c] << ": ";
        std::cout << constraintsPtrVec_[c]->m << " [";
        std::cout << iter << " ";
        iter += constraintsPtrVec_[c]->m;
        std::cout << iter << "]\n";
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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of init in get_starting_point!");
    }

    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_starting_point!");
    }

    if (x0.size() != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error x0.size() != numVars in get_starting_point!");
    }

    for ( Index i = 0; i < n; i++ ) {
        x[i] = x0(i);
    }

    return true;
}
// [TNLP_get_starting_point]

// [TNLP_update_minimal_cost_solution]
bool Optimizer::update_minimal_cost_solution(
    Index         n,
    const Number* x,
    Number        obj_value
) 
{
    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in update_minimal_cost_solution!");
    }

    if (obj_value < currentMinimalCost) {
        minimalCostSolution.resize(n);
        for ( Index i = 0; i < n; i++ ) {
            minimalCostSolution(i) = x[i];
        }
        currentMinimalCost = obj_value;
    }

    return true;
}
// [TNLP_update_minimal_cost_solution]

// [eval_hess_f]
// returns the hessian of the objective
bool Optimizer::eval_hess_f(
    Index         n,
    const Number* x,
    bool          new_x,
    MatX&         hess_f
) 
{
    throw std::invalid_argument("Objective function hessian is not implemented! Disable hessian please.");
    return false;
}
// [eval_hess_f]

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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_g!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in eval_g!");
    }

    // fill in a Eigen Vector instance of decision variables
    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    Index iter = 0;
    bool ifFeasibleCurrIter = true;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // compute constraints
        if (output_computation_time) {
            start_time = std::chrono::high_resolution_clock::now();
        }
        
        try {
            constraintsPtrVec_[c]->compute(z, false);
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in eval_g: " + constraintsNameVec_[c] + "! Check previous error message.");
        }

        if (output_computation_time) {
            end_time = std::chrono::high_resolution_clock::now();
            std::cout << "eval_g compute time for " << constraintsNameVec_[c] 
                      << " is " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() 
                      << " microseconds.\n";
        }

        // test if constraints are feasible
        if ((constraintsPtrVec_[c]->g - constraintsPtrVec_[c]->g_lb).minCoeff() < 0 || 
            (constraintsPtrVec_[c]->g_ub - constraintsPtrVec_[c]->g).minCoeff() < 0) {
            ifFeasibleCurrIter = false;
        }

        // fill in constraints
        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g[iter] = constraintsPtrVec_[c]->g(i);
            iter++;
        }
    }

    if (ifFeasibleCurrIter) {
        lastFeasibleSolution = z;
    }

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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_jac_g!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in eval_jac_g!");
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

        Index iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            // compute constraints
            if (output_computation_time) {
                start_time = std::chrono::high_resolution_clock::now();
            }

            try {
                constraintsPtrVec_[c]->compute(z, true);
            }
            catch (std::exception& e) {
                std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
                std::cout << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_jac_g in: " + constraintsNameVec_[c] + "! Check previous error message.");
            }

            if (output_computation_time) {
                end_time = std::chrono::high_resolution_clock::now();
                std::cout << "eval_jac_g compute time for " << constraintsNameVec_[c] 
                          << " is " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() 
                          << " microseconds.\n";
            }

            // fill in constraints
            for ( Index i = 0; i < constraintsPtrVec_[c]->pg_pz.rows(); i++ ) {
                for ( Index j = 0; j < constraintsPtrVec_[c]->pg_pz.cols(); j++ ) {
                    values[iter] = constraintsPtrVec_[c]->pg_pz(i, j);
                    iter++;
                }
            }
        }
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
    if (!enable_hessian) {
        return false;
    }

    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_h!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in eval_h!");
    }
    if (nele_hess != n * (n + 1) / 2) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of nele_hess in eval_h!");
    }

    if (values == NULL) {
        /* return the structure */
        /* the hessian for this problem is actually dense */
        Index iter = 0;
        for(Index i = 0; i < n; i++){
            for(Index j = i; j < n; j++){
                iRow[iter] = i;
                jCol[iter] = j;
                iter++;
            }
        }
    }
    else {
        // fill in a Eigen Vector instance of decision variables
        VecX z(n);
        for ( Index i = 0; i < n; i++ ) {
            z(i) = x[i];
        }

        MatX hessian(n, n);
        hessian.setZero();

        MatX hess_f(n, n);
        hess_f.setZero();

        if (obj_factor != 0) {
            if (output_computation_time) {
                start_time = std::chrono::high_resolution_clock::now();
            }

            try {
                eval_hess_f(n, x, new_x, hess_f);
            }
            catch (std::exception& e) {
                std::cerr << "Error in eval_hess_f: " << std::endl;
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_hess_f! Check previous error message.");
            }

            if (output_computation_time) {
                end_time = std::chrono::high_resolution_clock::now();
                std::cout << "eval_hess_f compute time is" << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() 
                          << " microseconds.\n";
            }

            hessian += obj_factor * hess_f;
        }

        Index iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            // if lambda is zero for this constraint, skip (make CheckDerivative faster)
            bool lambdaAllZeroForThisConstraint = true;
            for (Index innerIter = iter; innerIter < iter + constraintsPtrVec_[c]->m; innerIter++) {
                if (lambda[innerIter] != 0) {
                    lambdaAllZeroForThisConstraint = false;
                    break;
                }
            }
            if (lambdaAllZeroForThisConstraint) {
                iter += constraintsPtrVec_[c]->m;
                continue;
            } 

            // compute constraints
            if (output_computation_time) {
                start_time = std::chrono::high_resolution_clock::now();
            }

            try {
                constraintsPtrVec_[c]->compute(z, true, true);
            }
            catch (std::exception& e) {
                std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
                std::cout << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_h in: " + constraintsNameVec_[c] + "! Check previous error message.");
            }

            if (output_computation_time) {
                end_time = std::chrono::high_resolution_clock::now();
                std::cout << "eval_hessian_g compute time for " << constraintsNameVec_[c] 
                          << " is " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() 
                          << " microseconds.\n";
            }

            if (constraintsPtrVec_[c]->pg_pz_pz.size() != constraintsPtrVec_[c]->m) {
                std::cerr << "Error in " 
                          << constraintsNameVec_[c] 
                          << ": " 
                          << constraintsPtrVec_[c]->pg_pz_pz.size() 
                          << " != "
                          << constraintsPtrVec_[c]->m << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ have wrong sizes!");
            }

            for (Index i = 0; i < constraintsPtrVec_[c]->m; i++) {
                if (constraintsPtrVec_[c]->pg_pz_pz(i).rows() != n || 
                    constraintsPtrVec_[c]->pg_pz_pz(i).cols() != n) {
                    std::cerr << "Error in " 
                              << constraintsNameVec_[c] 
                              << ": " 
                              << constraintsPtrVec_[c]->pg_pz_pz(i).rows() 
                              << "x" 
                              << constraintsPtrVec_[c]->pg_pz_pz(i).cols() 
                              << " != " 
                              << n 
                              << "x" 
                              << n << std::endl;
                    THROW_EXCEPTION(IpoptException, "*** Error constraintsPtrVec_ have wrong sizes!");
                }

                hessian += lambda[iter] * constraintsPtrVec_[c]->pg_pz_pz(i);
                iter++;
            }
        }

        iter = 0;
        for(Index i = 0; i < n; i++){
            for(Index j = i; j < n; j++){
                values[iter] = hessian(i, j);
                iter++;
            }
        }

        // auto end_time = std::chrono::high_resolution_clock::now();
        // std::cout << "eval_h time: " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << " microseconds.\n";
    }

    return true;
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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in finalize_solution!");
    }
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in finalize_solution!");
    }

    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.

    // store the solution
    solution.resize(n);
    for( Index i = 0; i < n; i++ ) {
        solution(i) = x[i];
    }

    // re-evaluate constraints to update values in each constraint instances
    g_copy.resize(m);

    if(currentMinimalCost < obj_value) {
        Number x_copy[n];
        for( Index i = 0; i < n; i++ ) {
            x_copy[i] = solution(i);
        }

        eval_g(n, x_copy, false, m, g_copy.data());
        summarize_constraints(m, g, false);

        if (ifFeasible) {
            std::cout << "The final solution returned by ipopt is not the minimal cost solution! " << std::endl;
            std::cout << "Switch to the minimal cost solution recorded before." << std::endl;
            solution = minimalCostSolution;
            obj_value_copy = currentMinimalCost;
        }
    }
    else {
        obj_value_copy = obj_value;
        eval_g(n, x, false, m, g_copy.data());
        summarize_constraints(m, g);
    }

    if ((!ifFeasible) && lastFeasibleSolution.size() == n) {
        ifFeasible = true;
        solution = lastFeasibleSolution;
        std::cout << "Solution is not feasible but we have found one feasible solution before. " << std::endl;
        std::cout << "Switch to the last feasible solution recorded before." << std::endl;

        Number x_new[n];
        for( Index i = 0; i < n; i++ ) {
            x_new[i] = solution(i);
        }

        // re-evaluate constraints to update values in each constraint instances
        eval_f(n, x_new, false, obj_value);
        eval_g(n, x_new, false, m, g_copy.data());
        obj_value_copy = obj_value;
    }

    std::cout << "Objective value: " << obj_value_copy << std::endl;
}
// [TNLP_finalize_solution]

// [TNLP_summarize_constraints]
void Optimizer::summarize_constraints(
    Index                      m,
    const Number*              g,
    const bool                 verbose
) 
{
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in summarize_constraints!");
    }

    if (verbose) std::cout << "Constraint violation report:" << std::endl;

    Index iter = 0;
    ifFeasible = true;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // find where the maximum constraint violation is
        Number max_constr_violation = 0;
        Index max_constr_violation_id1 = 0;
        Index max_constr_violation_id2 = 0;
        for (Index i = 0; i < constraintsPtrVec_[c]->m; i++) {
            Number constr_violation = std::max(constraintsPtrVec_[c]->g_lb(i) - g[iter], 
                                               g[iter] - constraintsPtrVec_[c]->g_ub(i));
           
            if (constr_violation > max_constr_violation) {
                max_constr_violation_id1 = i;
                max_constr_violation_id2 = iter;
                max_constr_violation = constr_violation;
            }

            if (constr_violation > final_constr_violation) {
                final_constr_violation = constr_violation;
            }

            iter++;
        }

        // report constraint violation
        if (max_constr_violation > 0) {
            if (verbose) {
                std::cout << constraintsNameVec_[c] << ": " << max_constr_violation << std::endl;
                std::cout << "    range: [" << constraintsPtrVec_[c]->g_lb[max_constr_violation_id1] 
                                            << ", " 
                                            << constraintsPtrVec_[c]->g_ub[max_constr_violation_id1] 
                          << "], value: "   << g[max_constr_violation_id2] << std::endl;
            }

            constraintsPtrVec_[c]->print_violation_info();
                    
            ifFeasible = false;
        }
    }

    if (verbose) std::cout << "Total constraint violation: " << final_constr_violation << std::endl;
}
// [TNLP_summarize_constraints]

}; // namespace IDTO