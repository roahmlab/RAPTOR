#include "DualArmourOptimizer.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

bool DualArmourOptimizer::set_parameters(
    const VecX& q_des_input,
    Number t_plan_input,
    const std::shared_ptr<RobotInfo>& robotInfoPtr1_input,
    const std::shared_ptr<BezierCurveInterval>& trajPtr1_input,
    const std::shared_ptr<PZDynamics>& dynPtr1_input,
    const std::shared_ptr<RobotInfo>& robotInfoPtr2_input,
    const std::shared_ptr<BezierCurveInterval>& trajPtr2_input,
    const std::shared_ptr<PZDynamics>& dynPtr2_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input
 ) 
 {
    enable_hessian = false;

    if (q_des_input.size() != 2 * NUM_FACTORS) {
        throw std::invalid_argument("q_des_input.size() != 2 * NUM_FACTORS");
    }

    armourOptPtr1 = std::make_shared<ArmourOptimizer>();
    armourOptPtr2 = std::make_shared<ArmourOptimizer>();

    armourOptPtr1->set_parameters(q_des_input.head(NUM_FACTORS), 
                                  t_plan_input, 
                                  robotInfoPtr1_input, 
                                  trajPtr1_input, 
                                  dynPtr1_input, 
                                  boxCenters_input, 
                                  boxOrientation_input, 
                                  boxSize_input);
    armourOptPtr2->set_parameters(q_des_input.tail(NUM_FACTORS),
                                  t_plan_input,
                                  robotInfoPtr2_input,
                                  trajPtr2_input,
                                  dynPtr2_input,
                                  boxCenters_input,
                                  boxOrientation_input,
                                  boxSize_input);

    return true;
}

bool DualArmourOptimizer::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    numCons = 0;

    Index m1 = 0;
    armourOptPtr1->get_nlp_info(n, m1, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m1;

    Index m2 = 0;
    armourOptPtr2->get_nlp_info(n, m2, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m2;

    // TODO: numCons += arm-arm collision constraints?
    
    m = numCons;

    // The problem described 2 * NUM_FACTORS variables, since there are two arms
    numVars = 2 * NUM_FACTORS;
    n = 2 * NUM_FACTORS;

    // the nonzero structure of the Jacobian is the same for both arms, so they are stored separately
    nnz_jac_g = 2 * m * NUM_FACTORS; // TODO: + nnz for arm-arm collision constraint jacobians?

    nnz_h_lag = n * n;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool DualArmourOptimizer::get_bounds_info(
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
    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_bounds_info!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in get_bounds_info!");
    }

    armourOptPtr1->get_bounds_info(NUM_FACTORS, 
                                   x_l, 
                                   x_u, 
                                   armourOptPtr1->numCons, 
                                   g_l, 
                                   g_u);
    armourOptPtr2->get_bounds_info(NUM_FACTORS, 
                                   x_l + NUM_FACTORS, 
                                   x_u + NUM_FACTORS, 
                                   armourOptPtr2->numCons, 
                                   g_l + armourOptPtr1->numCons, 
                                   g_u + armourOptPtr1->numCons);

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool DualArmourOptimizer::get_starting_point(
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
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of init in get_starting_point!");
    }

    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_starting_point!");
    }

    for( Index i = 0; i < n; i++ ) {
        // initialize to zero
        x[i] = 0.0;

        // try to avoid local minimum
        // x[i] = min(max((q_des[i] - trajPtr_->q0[i]) / k_range[i], -0.5), 0.5);
    }

    return true;
}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool DualArmourOptimizer::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_f!");
    }

    Number obj_value1 = 0.0;
    armourOptPtr1->eval_f(NUM_FACTORS, x, new_x, obj_value1);

    Number obj_value2 = 0.0;
    armourOptPtr2->eval_f(NUM_FACTORS, x + NUM_FACTORS, new_x, obj_value2);

    obj_value = obj_value1 + obj_value2;

    update_minimal_cost_solution(n, Utils::initializeEigenVectorFromArray(x, n), new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DualArmourOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_grad_f!");
    }

    armourOptPtr1->eval_grad_f(NUM_FACTORS, x, new_x, grad_f);
    armourOptPtr2->eval_grad_f(NUM_FACTORS, x + NUM_FACTORS, new_x, grad_f + NUM_FACTORS);

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool DualArmourOptimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_g!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "Error wrong value of m in eval_g!");
    }

    try {
        armourOptPtr1->eval_g(NUM_FACTORS, x, new_x, armourOptPtr1->numCons, g);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the first robot!");
    }

    try {
        armourOptPtr2->eval_g(NUM_FACTORS, x + NUM_FACTORS, new_x, armourOptPtr2->numCons, g + armourOptPtr1->numCons);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the second robot!");
    }

    // TODO: update g for arm-arm collision constraints, note that sphere locations have been computed in eval_g above
    // start location: g + armourOptPtr1->numCons + armourOptPtr2->numCons

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool DualArmourOptimizer::eval_jac_g(
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
    if(n != 2 * NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_jac_g!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "Error wrong value of m in eval_jac_g!");
    }
        
    if( values == NULL ) {
        // return the structure of the Jacobian
        // this particular Jacobian is sparse

        Index idx = 0;

        // first robot
        for(Index i = 0; i < armourOptPtr1->numCons; i++){
            for(Index j = 0; j < NUM_FACTORS; j++){
                iRow[idx] = i;
                jCol[idx] = j;
                idx++;
            }
        }

        // second robot
        for (Index i = 0; i < armourOptPtr2->numCons; i++) {
            for (Index j = 0; j < NUM_FACTORS; j++) {
                iRow[idx] = i + armourOptPtr1->numCons;
                jCol[idx] = j + NUM_FACTORS;
                idx++;
            }
        }

        // TODO: update iRow and jCol for arm-arm collision constraints
    }
    else {
        try {
            armourOptPtr1->eval_jac_g(NUM_FACTORS, x, new_x, armourOptPtr1->numCons, nele_jac, iRow, jCol, values);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the first robot!");
        }

        try {
            armourOptPtr2->eval_jac_g(NUM_FACTORS, x + NUM_FACTORS, new_x, armourOptPtr2->numCons, nele_jac, iRow, jCol, values + armourOptPtr1->numCons * NUM_FACTORS);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the second robot!");
        }

        // TODO: update values for arm-arm collision constraints jacobians
    }

    return true;
}
// [TNLP_eval_jac_g]

void DualArmourOptimizer::summarize_constraints(
    Index                      m,
    const Number*              g,
    const bool                 verbose
) 
{
    ifFeasible = true;
    final_constr_violation = 0;

    // first robot
    armourOptPtr1->summarize_constraints(armourOptPtr1->numCons, g, verbose);

    ifFeasible &= armourOptPtr1->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, armourOptPtr1->final_constr_violation);

    // second robot
    armourOptPtr2->summarize_constraints(armourOptPtr2->numCons, g + armourOptPtr1->numCons, verbose);

    ifFeasible &= armourOptPtr2->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, armourOptPtr2->final_constr_violation);

    // TODO: update ifFeasible and final_constr_violation for arm-arm collision constraints
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR