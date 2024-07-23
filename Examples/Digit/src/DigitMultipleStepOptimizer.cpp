#include "DigitMultipleStepOptimizer.h"

namespace IDTO {
namespace Digit {

// // constructor
// DigitMultipleStepOptimizer::DigitMultipleStepOptimizer()
// {
// }


// // destructor
// DigitMultipleStepOptimizer::~DigitMultipleStepOptimizer()
// {
// }

bool DigitMultipleStepOptimizer::set_parameters(
    const int NSteps_input,
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const TimeDiscretization time_discretization_input,
    const int degree_input,
    const Model& model_input, 
    const Eigen::VectorXi& jtype_input,
    const std::vector<GaitParameters>& gps_input
) 
{
    if (gps_input.size() != NSteps_input) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong size of gps_input in set_parameters!");
    }

    x0 = x0_input;
    
    stepOptVec_.reserve(NSteps_input);
    for (int i = 0; i < NSteps_input; i++) {
        stepOptVec_.push_back(std::make_shared<DigitSingleStepOptimizer>());

        char stanceLeg = (i % 2 == 0) ? 
            'L' : 
            'R';
        Transform stance_foot_T_des = (i % 2 == 0) ? 
            Transform(3, -M_PI / 2) : 
            Transform(3, M_PI / 2);

        stepOptVec_[i]->set_parameters(x0_input,
                                       T_input,
                                       N_input,
                                       time_discretization_input,
                                       degree_input,
                                       model_input, 
                                       jtype_input,
                                       gps_input[i],
                                       stanceLeg,
                                       stance_foot_T_des,
                                       false);
    }

    const frictionParams FRICTION_PARAMS(MU, GAMMA, FOOT_WIDTH, FOOT_LENGTH);

    periodConsVec_.reserve(NSteps_input);
    for (int i = 0; i < NSteps_input - 1; i++) {
        periodConsVec_.push_back(
            std::make_shared<DigitMultipleStepPeriodicityConstraints>(
                stepOptVec_[i]->trajPtr_,
                stepOptVec_[i + 1]->trajPtr_,
                stepOptVec_[i]->dcidPtr_,
                stepOptVec_[i + 1]->dcidPtr_,
                FRICTION_PARAMS
            ));
    }

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
bool DigitMultipleStepOptimizer::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    numVars = 0;
    numCons = 0;

    n_local.resize(stepOptVec_.size());
    n_position.resize(stepOptVec_.size() + 1);
    m_local.resize(stepOptVec_.size());
    m_position.resize(stepOptVec_.size() + 1);

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->get_nlp_info(n_local[i], m_local[i], nnz_jac_g, nnz_h_lag, index_style);

        n_position[i] = numVars;
        m_position[i] = numCons;

        numVars += n_local[i];
        numCons += m_local[i];

        n_position[i + 1] = numVars;
        m_position[i + 1] = numCons;
    }

    for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) {
        m_position[stepOptVec_.size() + i] = numCons;

        numCons += periodConsVec_[i]->m;

        m_position[stepOptVec_.size() + i + 1] = numCons;
    }

    n = numVars;
    m = numCons;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool DigitMultipleStepOptimizer::get_bounds_info(
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

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        std::cout << "gait " << i << " bounds infomation:" << std::endl;
        stepOptVec_[i]->get_bounds_info(n_local[i], x_l + n_position[i], x_u + n_position[i], 
                                        m_local[i], g_l + m_position[i], g_u + m_position[i]);
    }

    for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) {
        const Index start_pos = m_position[stepOptVec_.size() + i];
        const Index end_pos = m_position[stepOptVec_.size() + i + 1];

        std::cout << "gait " << i << " - " << i + 1 << " continuous constraint: "
                  << periodConsVec_[i]->m 
                  << " [" << start_pos << " " << end_pos << "]" << std::endl;

        periodConsVec_[i]->compute_bounds();

        for ( Index j = start_pos; j < end_pos; j++ ) {
            g_l[j] = periodConsVec_[i]->g_lb(j - start_pos);
            g_u[j] = periodConsVec_[i]->g_ub(j - start_pos);
        }
    }

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_eval_f]
// returns the value of the objective function
bool DigitMultipleStepOptimizer::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    obj_value = 0;

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        Number obj_value_local;
        stepOptVec_[i]->eval_f(n_local[i], x + n_position[i], new_x, obj_value_local);

        obj_value += obj_value_local;
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DigitMultipleStepOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->eval_grad_f(n_local[i], x + n_position[i], new_x, grad_f + n_position[i]);
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool DigitMultipleStepOptimizer::eval_g(
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

    ifFeasibleCurrIter = false;
    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->eval_g(n_local[i], x + n_position[i], new_x, m_local[i], g + m_position[i]);
        ifFeasibleCurrIter = ifFeasibleCurrIter & stepOptVec_[i]->ifFeasibleCurrIter;
    }

    for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) {
        VecX z_curr = Utils::initializeEigenVectorFromArray(x + n_position[i], n_local[i]);
        periodConsVec_[i]->compute(z_curr, false);

        const Index start_pos = m_position[stepOptVec_.size() + i];
        for ( Index j = 0; j < periodConsVec_[i]->m; j++ ) {
            g[start_pos + j] = periodConsVec_[i]->g(j);

            if (g[start_pos + j] < periodConsVec_[i]->g_lb(j) - constr_viol_tol || 
                g[start_pos + j] > periodConsVec_[i]->g_ub(j) + constr_viol_tol) {
                ifFeasibleCurrIter = false;
            }
        }
    }

    // the following code is for updating the optimal solution, directly copied from Optimizer.cpp
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
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

    return true;
}
// [TNLP_eval_g]

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool DigitMultipleStepOptimizer::eval_jac_g(
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
        // this particular Jacobian is dense in blocks
        Index idx = 0;
        for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
            for ( Index j = 0; j < m_local[i]; j++ ) {
                for ( Index k = 0; k < n_local[i]; k++ ) {
                    iRow[idx] = m_position[i] + j;
                    jCol[idx] = n_position[i] + k;
                    idx++;
                }
            }
        }

        for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) { 
            const Index start_pos = m_position[stepOptVec_.size() + i];
            const Index end_pos = m_position[stepOptVec_.size() + i + 1];

            for ( Index j = 0; j < periodConsVec_[i]->m; j++ ) {
                for ( Index k = 0; k < n_local[i]; k++ ) {
                    iRow[idx] = start_pos + j;
                    jCol[idx] = n_position[i] + k;
                    idx++;
                }
            }

            for ( Index j = 0; j < NUM_INDEPENDENT_JOINTS; j++ ) {
                for ( Index k = 0; k < n_local[i + 1]; k++ ) {
                    iRow[idx] = start_pos + NUM_JOINTS + NUM_DEPENDENT_JOINTS + j;
                    jCol[idx] = n_position[i + 1] + k;
                    idx++;
                }
            }

            for ( Index j = 0; j < NUM_INDEPENDENT_JOINTS; j++ ) {
                for ( Index k = 0; k < n_local[i + 1]; k++ ) {
                    iRow[idx] = start_pos + NUM_JOINTS + NUM_DEPENDENT_JOINTS + NUM_INDEPENDENT_JOINTS + j;
                    jCol[idx] = n_position[i + 1] + k;
                    idx++;
                }
            }
        }
    }
    else {
        Index idx = 0;
        for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
            stepOptVec_[i]->eval_jac_g(n_local[i], 
                                       x + n_position[i], 
                                       new_x, 
                                       m_local[i], 
                                       n_local[i] * m_local[i], 
                                       nullptr, 
                                       nullptr,
                                       values + idx);

            idx += n_local[i] * m_local[i];
        }

        for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) { 
            VecX z_curr = Utils::initializeEigenVectorFromArray(x + n_position[i], n_local[i]);
            periodConsVec_[i]->compute(z_curr, true);

            for ( Index j = 0; j < periodConsVec_[i]->m; j++ ) {
                for ( Index k = 0; k < n_local[i]; k++ ) {
                    values[idx] = periodConsVec_[i]->pg_pz(j, k);
                    idx++;
                }
            }

            for ( Index j = 0; j < NUM_INDEPENDENT_JOINTS; j++ ) {
                for ( Index k = 0; k < n_local[i + 1]; k++ ) {
                    values[idx] = periodConsVec_[i]->pg3_pz2(j, k);
                    idx++;
                }
            }

            for ( Index j = 0; j < NUM_INDEPENDENT_JOINTS; j++ ) {
                for ( Index k = 0; k < n_local[i + 1]; k++ ) {
                    values[idx] = periodConsVec_[i]->pg4_pz2(j, k);
                    idx++;
                }
            }
        }
    }

    return true;
}
// [TNLP_eval_jac_g]

// [TNLP_summarize_constraints]
void DigitMultipleStepOptimizer::summarize_constraints(
    Index                      m,
    const Number*              g,
    const bool                 verbose
) 
{
    if (m != numCons) {
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in summarize_constraints!");
    }

    if (verbose) std::cout << "Constraint violation report:" << std::endl;

    ifFeasible = true;
    final_constr_violation = 0;

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->summarize_constraints(m_local[i], g + m_position[i], verbose);

        ifFeasible = ifFeasible & stepOptVec_[i]->ifFeasible;
        final_constr_violation = std::max(final_constr_violation, stepOptVec_[i]->final_constr_violation);
    }

    for ( Index i = 0; i < stepOptVec_.size() - 1; i++ ) {
        VecX z_curr = solution.segment(n_position[i], n_local[i]);
        periodConsVec_[i]->compute(z_curr, false);

        Number max_constr_violation = 0;
        Index max_constr_violation_id = 0;
        const Index start_pos = m_position[stepOptVec_.size() + i];
        for ( Index j = 0; j < periodConsVec_[i]->m; j++ ) {
            Number constr_violation = std::max(periodConsVec_[i]->g_lb(j) - periodConsVec_[i]->g(j), 
                                               periodConsVec_[i]->g(j) - periodConsVec_[i]->g_ub(j));

            if (constr_violation > max_constr_violation) {
                max_constr_violation_id = j;
                max_constr_violation = constr_violation;
            }

            if (constr_violation > final_constr_violation) {
                final_constr_violation = constr_violation;
            }
        }

        if (max_constr_violation > constr_viol_tol) {
            if (verbose) {
                std::cout << "gait " << i << " - " << i + 1 << " continuous constraint: " 
                          << max_constr_violation << std::endl;
                std::cout << "    range: [" << periodConsVec_[i]->g_lb[max_constr_violation_id] 
                                            << ", " 
                                            << periodConsVec_[i]->g_ub[max_constr_violation_id] 
                          << "], value: "   << periodConsVec_[i]->g(max_constr_violation_id) << std::endl;
            }

            periodConsVec_[i]->print_violation_info();
                    
            ifFeasible = false;
        }
    }

    if (verbose) std::cout << "Total constraint violation: " << final_constr_violation << std::endl;
}
// [TNLP_summarize_constraints]

}; // namespace Digit
}; // namespace IDTO