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

        char stanceLeg = (i % 2 == 0) ? 'L' : 'R';

        stepOptVec_[i]->set_parameters(x0_input,
                                       T_input,
                                       N_input,
                                       time_discretization_input,
                                       degree_input,
                                       model_input, 
                                       jtype_input,
                                       gps_input[i],
                                       stanceLeg,
                                       Transform(3, -M_PI / 2),
                                       false);
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
    m_local.resize(stepOptVec_.size());

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->get_nlp_info(n_local[i], m_local[i], nnz_jac_g, nnz_h_lag, index_style);
        numVars += n_local[i];
        numCons += m_local[i];
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

    Index n_offset = 0;
    Index m_offset = 0;

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        std::cout << "gait " << i << " bounds infomation:" << std::endl;
        stepOptVec_[i]->get_bounds_info(n_local[i], x_l + n_offset, x_u + n_offset, 
                                        m_local[i], g_l + m_offset, g_u + m_offset);

        n_offset += n_local[i]; 
        m_offset += m_local[i];
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
    Index n_offset = 0;

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        Number obj_value_local;
        stepOptVec_[i]->eval_f(n_local[i], x + n_offset, new_x, obj_value_local);

        obj_value += obj_value_local;
        n_offset += n_local[i];
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

    Index n_offset = 0;

    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->eval_grad_f(n_local[i], x + n_offset, new_x, grad_f + n_offset);
        n_offset += n_local[i];
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

    Index n_offset = 0;
    Index m_offset = 0;

    ifFeasibleCurrIter = false;
    for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
        stepOptVec_[i]->eval_g(n_local[i], x + n_offset, new_x, m_local[i], g + m_offset);
        ifFeasibleCurrIter = ifFeasibleCurrIter & stepOptVec_[i]->ifFeasibleCurrIter;

        n_offset += n_local[i];
        m_offset += m_local[i];
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
        Index n_offset = 0;
        Index m_offset = 0;
        Index idx = 0;
        for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
            for ( Index j = 0; j < m_local[i]; j++ ) {
                for ( Index k = 0; k < n_local[i]; k++ ) {
                    iRow[idx] = m_offset + j;
                    jCol[idx] = n_offset + k;
                    idx++;
                }
            }
            n_offset += n_local[i];
            m_offset += m_local[i];
        }
    }
    else {
        Index n_offset = 0;
        Index m_offset = 0;
        Index v_offset = 0;

        for ( Index i = 0; i < stepOptVec_.size(); i++ ) {
            stepOptVec_[i]->eval_jac_g(n_local[i], 
                                       x + n_offset, 
                                       new_x, 
                                       m_local[i], 
                                       n_local[i] * m_local[i], 
                                       nullptr, 
                                       nullptr,
                                       values + v_offset);

            n_offset += n_local[i];
            m_offset += m_local[i];
            v_offset += n_local[i] * m_local[i];
        }
    }

    return true;
}
// [TNLP_eval_jac_g]

}; // namespace Digit
}; // namespace IDTO