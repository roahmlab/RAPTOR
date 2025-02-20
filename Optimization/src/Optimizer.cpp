#include "Optimizer.h"

namespace RAPTOR {

// [TNLP_reset]
// Method to reset the optimizer
void Optimizer::reset() {
    numVars = 0;
    numCons = 0;

    constraintsPtrVec_.clear();
    constraintsNameVec_.clear();

    ifFeasibleCurrIter = false;
    currentIpoptSolution = VecX();
    currentIpoptObjValue = std::numeric_limits<Number>::max();
    ifCurrentIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;
    optimalIpoptSolution = VecX();
    optimalIpoptObjValue = std::numeric_limits<Number>::max();
    ifOptimalIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;

    nlp_f_time = 0;
    nlp_f_count = 0;
    nlp_g_time = 0;
    nlp_g_count = 0;
    nlp_grad_f_time = 0;
    nlp_grad_f_count = 0;
    nlp_jac_g_time = 0;
    nlp_jac_g_count = 0;
    nlp_hess_l_time = 0;
    nlp_hess_l_count = 0;
}

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

    if (costsPtrVec_.size() != costsWeightVec_.size() ||
        costsPtrVec_.size() != costsNameVec_.size()) {
        THROW_EXCEPTION(IpoptException, "*** Error costsPtrVec_, costsWeightVec_, and costsNameVec_ have different sizes!");
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
    if (display_info && constraintsPtrVec_.size() > 0) {
        std::cout << "Dimension of each constraints and their locations: \n";
        iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            std::cout << constraintsNameVec_[c] << ": ";
            std::cout << constraintsPtrVec_[c]->m << " [";
            std::cout << iter << " ";
            iter += constraintsPtrVec_[c]->m;
            std::cout << iter << "]\n";
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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of init in get_starting_point!");
    }

    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_starting_point!");
    }

    if (x0.size() != numVars) {
        std::cerr << "x0.size() = " << x0.size() << std::endl;
        std::cerr << "numVars = " << numVars << std::endl;
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
    const VecX&   z,
    bool          new_x,
    Number        obj_value
) 
{
    if (n != numVars) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in update_minimal_cost_solution!");
    }

    // update status of the current solution
    if (new_x) { // directly assign currentIpoptSolution if this x has never been evaluated before
        currentIpoptSolution = z;
        currentIpoptObjValue = obj_value;
        ifCurrentIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;
    }
    else { // update currentIpoptSolution
        if (Utils::ifTwoVectorEqual(currentIpoptSolution, z, 0)) {
            if (ifCurrentIpoptFeasible == OptimizerConstants::FeasibleState::UNINITIALIZED) {
                THROW_EXCEPTION(IpoptException, "*** Error ifCurrentIpoptFeasible is not initialized!");
            }
            else { // this has been evaluated in eval_g, just need to update the cost
                currentIpoptObjValue = obj_value;
            }
        }
        else {
            currentIpoptSolution = z;
            currentIpoptObjValue = obj_value;
            ifCurrentIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;
        }
    }

    // update the status of the optimal solution
    if (ifCurrentIpoptFeasible == OptimizerConstants::FeasibleState::FEASIBLE &&
        currentIpoptObjValue < optimalIpoptObjValue) {
        optimalIpoptSolution = currentIpoptSolution;
        optimalIpoptObjValue = currentIpoptObjValue;
        ifOptimalIpoptFeasible = ifCurrentIpoptFeasible;
    }

    return true;
}
// [TNLP_update_minimal_cost_solution]

// [eval_f]
// returns the objective value
bool Optimizer::eval_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number&       obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    if (costsPtrVec_.size() == 0) {
        std::cerr << "You should add at least one cost function or implement your own eval_f!" << std::endl;
        THROW_EXCEPTION(IpoptException, "*** Error costsPtrVec_ is empty!");
    }

    start_time = std::chrono::high_resolution_clock::now();

    // fill in a Eigen Vector instance of decision variables
    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    obj_value = 0;
    for (Index f = 0; f < costsPtrVec_.size(); f++) {
        // compute cost functions
        try {
            costsPtrVec_[f]->compute(z, false);
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << costsNameVec_[f] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in eval_f: " + costsNameVec_[f] + "! Check previous error message.");
        }

        obj_value += costsWeightVec_[f] * costsPtrVec_[f]->f;
    }

    update_minimal_cost_solution(n, z, new_x, obj_value);

    end_time = std::chrono::high_resolution_clock::now();
    nlp_f_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    nlp_f_count++;

    return true;
}
// [eval_f]

// [eval_grad_f]
// returns the gradient of the objective
bool Optimizer::eval_grad_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number*       grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    if (costsPtrVec_.size() == 0) {
        std::cerr << "You should add at least one cost function or implement your own eval_f!" << std::endl;
        THROW_EXCEPTION(IpoptException, "*** Error costsPtrVec_ is empty!");
    }

    start_time = std::chrono::high_resolution_clock::now();

    // fill in a Eigen Vector instance of decision variables
    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    for (Index i = 0; i < n; i++) {
        grad_f[i] = 0;
    }

    for (Index f = 0; f < costsPtrVec_.size(); f++) {
        // compute cost functions
        try {
            costsPtrVec_[f]->compute(z, true);
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << costsNameVec_[f] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in eval_grad_f: " + costsNameVec_[f] + "! Check previous error message.");
        }

        for (Index i = 0; i < n; i++) {
            grad_f[i] += costsWeightVec_[f] * costsPtrVec_[f]->grad_f(i);
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    nlp_grad_f_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    nlp_grad_f_count++;

    return true;
}
// [eval_grad_f]

// [eval_hess_f]
// returns the hessian of the objective
bool Optimizer::eval_hess_f(
    Index         n,
    const Number* x,
    bool          new_x,
    MatX&         hess_f
) 
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    if (costsPtrVec_.size() == 0) {
        std::cerr << "You should add at least one cost function or implement your own eval_f!" << std::endl;
        THROW_EXCEPTION(IpoptException, "*** Error costsPtrVec_ is empty!");
    }

    hess_f = MatX::Zero(n, n);

    // fill in a Eigen Vector instance of decision variables
    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    for (Index f = 0; f < costsPtrVec_.size(); f++) {
        // compute cost functions
        try {
            costsPtrVec_[f]->compute(z, true, true);
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << costsNameVec_[f] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in eval_hess_f: " + costsNameVec_[f] + "! Check previous error message.");
        }

        hess_f += costsWeightVec_[f] * costsPtrVec_[f]->hess_f;
    }

    return true;
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

    start_time = std::chrono::high_resolution_clock::now();

    // fill in a Eigen Vector instance of decision variables
    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    Index iter = 0;
    ifFeasibleCurrIter = true;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // compute constraints
        try {
            constraintsPtrVec_[c]->compute(z, false);
        }
        catch (std::exception& e) {
            std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "*** Error in eval_g: " + constraintsNameVec_[c] + "! Check previous error message.");
        }

        // test if constraints are feasible
        if ((constraintsPtrVec_[c]->g - constraintsPtrVec_[c]->g_lb).minCoeff() < -constr_viol_tol || 
            (constraintsPtrVec_[c]->g_ub - constraintsPtrVec_[c]->g).minCoeff() < -constr_viol_tol) {
            ifFeasibleCurrIter = false;
        }

        // fill in constraints
        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g[iter] = constraintsPtrVec_[c]->g(i);
            iter++;
        }
    }

    // update status of the current solution
    if (new_x) { // directly assign currentIpoptSolution if this x has never been evaluated before
        currentIpoptSolution = z;
        currentIpoptObjValue = std::numeric_limits<Number>::max();
        ifCurrentIpoptFeasible = ifFeasibleCurrIter ? 
                                     OptimizerConstants::FeasibleState::FEASIBLE : 
                                     OptimizerConstants::FeasibleState::INFEASIBLE;
    }
    else { // update currentIpoptSolution
        if (Utils::ifTwoVectorEqual(currentIpoptSolution, z, 0)) {
            if (currentIpoptObjValue == std::numeric_limits<Number>::max()) {
                THROW_EXCEPTION(IpoptException, "*** Error currentIpoptObjValue is not initialized!");
            }
            else { // this has been evaluated in eval_f, just need to update the feasibility
                ifCurrentIpoptFeasible = ifFeasibleCurrIter ? 
                                             OptimizerConstants::FeasibleState::FEASIBLE : 
                                             OptimizerConstants::FeasibleState::INFEASIBLE;
            }
        }
        else {
            currentIpoptSolution = z;
            currentIpoptObjValue = std::numeric_limits<Number>::max();
            ifCurrentIpoptFeasible = ifFeasibleCurrIter ? 
                                         OptimizerConstants::FeasibleState::FEASIBLE : 
                                         OptimizerConstants::FeasibleState::INFEASIBLE;
        }
    }

    // update the status of the optimal solution
    if (ifCurrentIpoptFeasible == OptimizerConstants::FeasibleState::FEASIBLE &&
        currentIpoptObjValue < optimalIpoptObjValue) {
        optimalIpoptSolution = currentIpoptSolution;
        optimalIpoptObjValue = currentIpoptObjValue;
        ifOptimalIpoptFeasible = ifCurrentIpoptFeasible;
    }

    end_time = std::chrono::high_resolution_clock::now();
    nlp_g_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    nlp_g_count++;

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
        start_time = std::chrono::high_resolution_clock::now();

        // fill in a Eigen Vector instance of decision variables
        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        Index iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            // compute constraints
            try {
                constraintsPtrVec_[c]->compute(z, true);
            }
            catch (std::exception& e) {
                std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_jac_g in: " + constraintsNameVec_[c] + "! Check previous error message.");
            }

            // fill in constraints
            for ( Index i = 0; i < constraintsPtrVec_[c]->pg_pz.rows(); i++ ) {
                for ( Index j = 0; j < constraintsPtrVec_[c]->pg_pz.cols(); j++ ) {
                    values[iter] = constraintsPtrVec_[c]->pg_pz(i, j);
                    iter++;
                }
            }
        }

        end_time = std::chrono::high_resolution_clock::now();
        nlp_jac_g_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        nlp_jac_g_count++;
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
        /* the hessian for this problem is actually dense 
           but we only require the upper triangle part since hessian should be symmetric*/
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
        start_time = std::chrono::high_resolution_clock::now();

        // fill in a Eigen Vector instance of decision variables
        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        MatX hessian(n, n);
        hessian.setZero();

        MatX hess_f(n, n);
        hess_f.setZero();

        if (obj_factor != 0) {
            try {
                eval_hess_f(n, x, new_x, hess_f);
            }
            catch (std::exception& e) {
                std::cerr << "Error in eval_hess_f: " << std::endl;
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_hess_f! Check previous error message.");
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
            try {
                constraintsPtrVec_[c]->compute(z, true, true);
            }
            catch (std::exception& e) {
                std::cerr << "Error in " << constraintsNameVec_[c] << ": " << std::endl;
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "*** Error in eval_h in: " + constraintsNameVec_[c] + "! Check previous error message.");
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

        end_time = std::chrono::high_resolution_clock::now();
        nlp_hess_l_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        nlp_hess_l_count++;
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
    eval_f(n, x, true, obj_value_copy);
    eval_g(n, x, true, m, g_copy.data());
    summarize_constraints(m, g, false);

    bool recordedOptimalSolutionAvailable = 
        optimalIpoptSolution.size() == n &&
        optimalIpoptObjValue < std::numeric_limits<Number>::max() &&
        ifOptimalIpoptFeasible == OptimizerConstants::FeasibleState::FEASIBLE;

    if (recordedOptimalSolutionAvailable &&
        (!ifFeasible || optimalIpoptObjValue < obj_value)) {
        ifFeasible = true;
        solution = optimalIpoptSolution;
        obj_value_copy = optimalIpoptObjValue;

        Number x_new[n];
        for( Index i = 0; i < n; i++ ) {
            x_new[i] = solution(i);
        }
        
        // re-evaluate constraints to update values in each constraint instances
        eval_f(n, x_new, true, obj_value_copy);
        eval_g(n, x_new, true, m, g_copy.data());
        summarize_constraints(m, g_copy.data());
        
        if (display_info) {
            std::cout << "Solution is not feasible or optimal but we have found one optimal feasible solution before" 
                      << " with cost: " << optimalIpoptObjValue
                      << " and constraint violation: " << final_constr_violation << std::endl;
        }
    }
    else {
        summarize_constraints(m, g);
    }

    if (display_info) std::cout << "Objective value: " << obj_value_copy << std::endl;

    // clear the previous information
    ifFeasibleCurrIter = false;
    currentIpoptSolution = VecX();
    currentIpoptObjValue = std::numeric_limits<Number>::max();
    ifCurrentIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;
    optimalIpoptSolution = VecX();
    optimalIpoptObjValue = std::numeric_limits<Number>::max();
    ifOptimalIpoptFeasible = OptimizerConstants::FeasibleState::UNINITIALIZED;

    if (display_info) {
        // report the time taken for each function
        printf("  optimizer  :     t_wall      (avg)    n_eval\n");
        printf("      nlp_f  | %8.2fms (%7.2fus) %8d\n", nlp_f_time / 1000.0, nlp_f_time / nlp_f_count, nlp_f_count);
        printf("      nlp_g  | %8.2fms (%7.2fus) %8d\n", nlp_g_time / 1000.0, nlp_g_time / nlp_g_count, nlp_g_count);
        printf(" nlp_grad_f  | %8.2fms (%7.2fus) %8d\n", nlp_grad_f_time / 1000.0, nlp_grad_f_time / nlp_grad_f_count, nlp_grad_f_count);
        printf("  nlp_jac_g  | %8.2fms (%7.2fus) %8d\n", nlp_jac_g_time / 1000.0, nlp_jac_g_time / nlp_jac_g_count, nlp_jac_g_count);
        if (enable_hessian) {
            printf("  nlp_hess_l | %8.2fms (%7.2fus) %8d\n", nlp_hess_l_time / 1000.0, nlp_hess_l_time / nlp_hess_l_count, nlp_hess_l_count);
        }
        double total_time = nlp_f_time + nlp_g_time + nlp_grad_f_time + nlp_jac_g_time + nlp_hess_l_time;
        printf("      total  | %8.2fms (%7.2fms) %8d\n", total_time / 1000.0, total_time / 1000.0, 1);
    }

    // clear the timing information
    nlp_f_time = 0;
    nlp_f_count = 0;
    nlp_g_time = 0;
    nlp_g_count = 0;
    nlp_grad_f_time = 0;
    nlp_grad_f_count = 0;
    nlp_jac_g_time = 0;
    nlp_jac_g_count = 0;
    nlp_hess_l_time = 0;
    nlp_hess_l_count = 0;
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

    if (display_info && verbose && constraintsPtrVec_.size() > 0) {
        std::cout << "Constraint violation report:" << std::endl;
    }

    Index iter = 0;
    ifFeasible = true;
    final_constr_violation = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // find where the maximum constraint violation is
        Number max_constr_violation = 0;
        Index max_constr_violation_id1 = 0;
        Index max_constr_violation_id2 = 0;
        for (Index i = 0; i < constraintsPtrVec_[c]->m; i++) {
            Number constr_violation = fmax(constraintsPtrVec_[c]->g_lb(i) - g[iter], 
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
        if (max_constr_violation > constr_viol_tol) {
            if (display_info && verbose) {
                std::cout << constraintsNameVec_[c] << ": " << max_constr_violation << std::endl;
                std::cout << "    range: [" << constraintsPtrVec_[c]->g_lb[max_constr_violation_id1] 
                                            << ", " 
                                            << constraintsPtrVec_[c]->g_ub[max_constr_violation_id1] 
                          << "], value: "   << g[max_constr_violation_id2] << std::endl;
            }

            if (verbose) constraintsPtrVec_[c]->print_violation_info();
                    
            ifFeasible = false;
        }
    }

    if (display_info && verbose && constraintsPtrVec_.size() > 0) {
        std::cout << "Total constraint violation: " << final_constr_violation << std::endl;
    }
}
// [TNLP_summarize_constraints]

}; // namespace RAPTOR