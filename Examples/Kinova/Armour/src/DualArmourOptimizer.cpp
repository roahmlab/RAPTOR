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

    robotInfoPtr1_ = robotInfoPtr1_input;
    robotInfoPtr2_ = robotInfoPtr2_input;

    armourOptPtr1_ = std::make_shared<ArmourOptimizer>();
    armourOptPtr2_ = std::make_shared<ArmourOptimizer>();

    armourOptPtr1_->set_parameters(q_des_input.head(NUM_FACTORS), 
                                   t_plan_input, 
                                   robotInfoPtr1_input, 
                                   trajPtr1_input, 
                                   dynPtr1_input, 
                                   boxCenters_input, 
                                   boxOrientation_input, 
                                   boxSize_input);
    armourOptPtr2_->set_parameters(q_des_input.tail(NUM_FACTORS),
                                   t_plan_input,
                                   robotInfoPtr2_input,
                                   trajPtr2_input,
                                   dynPtr2_input,
                                   boxCenters_input,
                                   boxOrientation_input,
                                   boxSize_input);

    if (armourOptPtr1_->num_time_steps != armourOptPtr2_->num_time_steps) {
        throw std::invalid_argument("armourOptPtr1_->num_time_steps != armourOptPtr2_->num_time_steps");
    }

    tccPtrs.resize(armourOptPtr1_->num_time_steps);
    for (size_t i = 0; i < armourOptPtr1_->num_time_steps; i++) {
        tccPtrs[i] = std::make_shared<TaperedCapsuleCollision<2 * NUM_FACTORS>>();
    }

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
    armourOptPtr1_->get_nlp_info(n, m1, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m1;

    Index m2 = 0;
    armourOptPtr2_->get_nlp_info(n, m2, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m2;

    // arm-arm collision constraints
    const Index m_arm_arm = armourOptPtr1_->num_time_steps * robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size();
    numCons += m_arm_arm;
    
    m = numCons;

    // The problem described 2 * NUM_FACTORS variables, since there are two arms
    numVars = 2 * NUM_FACTORS;
    n = numVars;

    // the nonzero structure of the Jacobian is the same for both arms, so they are stored separately
    nnz_jac_g = 
        (armourOptPtr1_->numCons + armourOptPtr2_->numCons) * NUM_FACTORS + 
        m_arm_arm * 2 * NUM_FACTORS;

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

    armourOptPtr1_->get_bounds_info(NUM_FACTORS, 
                                    x_l, 
                                    x_u, 
                                    armourOptPtr1_->numCons, 
                                    g_l, 
                                    g_u);
    armourOptPtr2_->get_bounds_info(NUM_FACTORS, 
                                    x_l + NUM_FACTORS, 
                                    x_u + NUM_FACTORS, 
                                    armourOptPtr2_->numCons, 
                                    g_l + armourOptPtr1_->numCons, 
                                    g_u + armourOptPtr1_->numCons);

    const Index offset = armourOptPtr1_->numCons + armourOptPtr2_->numCons;
    for (Index i = 0; i < armourOptPtr1_->num_time_steps * robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size(); i++) {
        g_l[offset + i] = 0.0;
        g_u[offset + i] = 1e19;
    }

    g_lb_copy = Utils::initializeEigenVectorFromArray(g_l, m);
    g_ub_copy = Utils::initializeEigenVectorFromArray(g_u, m);

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
    armourOptPtr1_->eval_f(NUM_FACTORS, x, new_x, obj_value1);

    Number obj_value2 = 0.0;
    armourOptPtr2_->eval_f(NUM_FACTORS, x + NUM_FACTORS, new_x, obj_value2);

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

    armourOptPtr1_->eval_grad_f(NUM_FACTORS, x, new_x, grad_f);
    armourOptPtr2_->eval_grad_f(NUM_FACTORS, x + NUM_FACTORS, new_x, grad_f + NUM_FACTORS);

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
        armourOptPtr1_->eval_g(NUM_FACTORS, x, new_x, armourOptPtr1_->numCons, g);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the first robot!");
    }

    try {
        armourOptPtr2_->eval_g(NUM_FACTORS, x + NUM_FACTORS, new_x, armourOptPtr2_->numCons, g + armourOptPtr1_->numCons);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the second robot!");
    }

    // arm-arm collision constraints
    const Index offset = armourOptPtr1_->numCons + armourOptPtr2_->numCons;
    const auto& dynPtr1_ = armourOptPtr1_->dynPtr_;
    const auto& dynPtr2_ = armourOptPtr2_->dynPtr_;
    const size_t num_arm_arm_collision = robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size();

    try {
        Index i = 0;
        #pragma omp parallel for shared(armourOptPtr1_, armourOptPtr2_, robotInfoPtr1_, robotInfoPtr2_, x, g) private(i) schedule(dynamic)
        for(i = 0; i < armourOptPtr1_->num_time_steps; i++){
            Index j = 0;
            for (Index j1 = 0; j1 < robotInfoPtr1_->tc_begin_and_end.size(); j1++) {
                for (Index j2 = 0; j2 < robotInfoPtr2_->tc_begin_and_end.size(); j2++) {
                    const size_t tc1_begin_index = robotInfoPtr1_->tc_begin_and_end[j1].first;
                    const size_t tc1_end_index   = robotInfoPtr1_->tc_begin_and_end[j1].second;
                    const size_t tc2_begin_index = robotInfoPtr2_->tc_begin_and_end[j2].first;
                    const size_t tc2_end_index   = robotInfoPtr2_->tc_begin_and_end[j2].second;

                    const Vec3& tc1_sphere_1 = armourOptPtr1_->sphere_locations(i, tc1_begin_index);
                    const Vec3& tc1_sphere_2 = armourOptPtr1_->sphere_locations(i, tc1_end_index);
                    const Vec3& tc2_sphere_1 = armourOptPtr2_->sphere_locations(i, tc2_begin_index);
                    const Vec3& tc2_sphere_2 = armourOptPtr2_->sphere_locations(i, tc2_end_index);

                    const double tc1_sphere_1_radius = dynPtr1_->sphere_radii(tc1_begin_index, i);
                    const double tc1_sphere_2_radius = dynPtr1_->sphere_radii(tc1_end_index, i);
                    const double tc2_sphere_1_radius = dynPtr2_->sphere_radii(tc2_begin_index, i);
                    const double tc2_sphere_2_radius = dynPtr2_->sphere_radii(tc2_end_index, i);

                    g[i * num_arm_arm_collision + j + offset] = tccPtrs[i]->computeDistance(
                        tc1_sphere_1, tc1_sphere_2, 
                        tc2_sphere_1, tc2_sphere_2,
                        tc1_sphere_1_radius, tc1_sphere_2_radius, 
                        tc2_sphere_1_radius, tc2_sphere_2_radius);
                    j++;
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the arm-arm collision constraints!");
    }

    // update status of the current solution 
    // originally computed in Optimizer.cpp but we have overwitten eval_g, so have to manually update it here
    ifFeasibleCurrIter = true;
    for (Index i = 0; i < m; i++) {
        if (std::isnan(g[i])) {
            std::cerr << "g[" << i << "] is nan!" << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_g!");
        }
        
        // test if constraints are feasible
        if (g[i] - g_lb_copy[i] < -constr_viol_tol || 
            g_ub_copy[i] - g[i] < -constr_viol_tol) {
            ifFeasibleCurrIter = false;
            break;
        }
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);
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
        for(Index i = 0; i < armourOptPtr1_->numCons; i++){
            for(Index j = 0; j < NUM_FACTORS; j++){
                iRow[idx] = i;
                jCol[idx] = j;
                idx++;
            }
        }

        // second robot
        for (Index i = 0; i < armourOptPtr2_->numCons; i++) {
            for (Index j = 0; j < NUM_FACTORS; j++) {
                iRow[idx] = i + armourOptPtr1_->numCons;
                jCol[idx] = j + NUM_FACTORS;
                idx++;
            }
        }

        // update iRow and jCol for arm-arm collision constraints
        Index offset = armourOptPtr1_->numCons + armourOptPtr2_->numCons;
        const size_t num_arm_arm_collision = robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size();
        for (Index i = 0; i < armourOptPtr1_->num_time_steps; i++) {
            Index j = 0;
            for (Index j1 = 0; j1 < robotInfoPtr1_->tc_begin_and_end.size(); j1++) {
                for (Index j2 = 0; j2 < robotInfoPtr2_->tc_begin_and_end.size(); j2++) {
                    for (Index k = 0; k < 2 * NUM_FACTORS; k++) {
                        iRow[idx] = i * num_arm_arm_collision + j + offset;
                        jCol[idx] = k;
                        idx++;
                    }
                    j++;
                }
            }
        }
    }
    else {
        try {
            armourOptPtr1_->eval_jac_g(NUM_FACTORS, x, new_x, armourOptPtr1_->numCons, nele_jac, iRow, jCol, values);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the first robot!");
        }

        try {
            armourOptPtr2_->eval_jac_g(NUM_FACTORS, x + NUM_FACTORS, new_x, armourOptPtr2_->numCons, nele_jac, iRow, jCol, values + armourOptPtr1_->numCons * NUM_FACTORS);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the second robot!");
        }

        // arm-arm collision constraints
        const Index offset = (armourOptPtr1_->numCons + armourOptPtr2_->numCons) * NUM_FACTORS;
        const auto& dynPtr1_ = armourOptPtr1_->dynPtr_;
        const auto& dynPtr2_ = armourOptPtr2_->dynPtr_;
        const size_t num_arm_arm_collision = robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size();

        try {
            Index i = 0;
            #pragma omp parallel for shared(armourOptPtr1_, armourOptPtr2_, robotInfoPtr1_, robotInfoPtr2_, x, values) private(i) schedule(dynamic)
            for(i = 0; i < armourOptPtr1_->num_time_steps; i++){
                Index j = 0;
                for (Index j1 = 0; j1 < robotInfoPtr1_->tc_begin_and_end.size(); j1++) {
                    for (Index j2 = 0; j2 < robotInfoPtr2_->tc_begin_and_end.size(); j2++) {
                        const size_t tc1_begin_index = robotInfoPtr1_->tc_begin_and_end[j1].first;
                        const size_t tc1_end_index = robotInfoPtr1_->tc_begin_and_end[j1].second;
                        const size_t tc2_begin_index = robotInfoPtr2_->tc_begin_and_end[j2].first;
                        const size_t tc2_end_index = robotInfoPtr2_->tc_begin_and_end[j2].second;

                        const Vec3& tc1_sphere_1 = armourOptPtr1_->sphere_locations(i, tc1_begin_index);
                        const Vec3& tc1_sphere_2 = armourOptPtr1_->sphere_locations(i, tc1_end_index);
                        const Vec3& tc2_sphere_1 = armourOptPtr2_->sphere_locations(i, tc2_begin_index);
                        const Vec3& tc2_sphere_2 = armourOptPtr2_->sphere_locations(i, tc2_end_index);

                        const double tc1_sphere_1_radius = dynPtr1_->sphere_radii(tc1_begin_index, i);
                        const double tc1_sphere_2_radius = dynPtr1_->sphere_radii(tc1_end_index, i);
                        const double tc2_sphere_1_radius = dynPtr2_->sphere_radii(tc2_begin_index, i);
                        const double tc2_sphere_2_radius = dynPtr2_->sphere_radii(tc2_end_index, i);

                        Eigen::Matrix<double, 3, 2 * NUM_FACTORS> dk_tc1_sphere_1;
                        Eigen::Matrix<double, 3, 2 * NUM_FACTORS> dk_tc1_sphere_2;
                        Eigen::Matrix<double, 3, 2 * NUM_FACTORS> dk_tc2_sphere_1;
                        Eigen::Matrix<double, 3, 2 * NUM_FACTORS> dk_tc2_sphere_2;
                        dk_tc1_sphere_1.setZero();
                        dk_tc1_sphere_2.setZero();
                        dk_tc2_sphere_1.setZero();
                        dk_tc2_sphere_2.setZero();

                        dk_tc1_sphere_1.leftCols(NUM_FACTORS) = armourOptPtr1_->sphere_gradient(i, tc1_begin_index);
                        dk_tc1_sphere_2.leftCols(NUM_FACTORS) = armourOptPtr1_->sphere_gradient(i, tc1_end_index);
                        dk_tc2_sphere_1.rightCols(NUM_FACTORS) = armourOptPtr2_->sphere_gradient(i, tc2_begin_index);
                        dk_tc2_sphere_2.rightCols(NUM_FACTORS) = armourOptPtr2_->sphere_gradient(i, tc2_end_index);

                        Eigen::Vector<double, 2 * NUM_FACTORS> dk_distances;
                        const double distance = tccPtrs[i]->computeDistance(
                            tc1_sphere_1, tc1_sphere_2, 
                            tc2_sphere_1, tc2_sphere_2,
                            dk_tc1_sphere_1, dk_tc1_sphere_2, 
                            dk_tc2_sphere_1, dk_tc2_sphere_2,
                            tc1_sphere_1_radius, tc1_sphere_2_radius, 
                            tc2_sphere_1_radius, tc2_sphere_2_radius,
                            dk_distances);

                        std::memcpy(values + (i * num_arm_arm_collision + j) * 2 * NUM_FACTORS + offset, dk_distances.data(), 2 * NUM_FACTORS * sizeof(Number));

                        j++;
                    }
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the arm-arm collision constraints!");
        }
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
    armourOptPtr1_->summarize_constraints(armourOptPtr1_->numCons, g, verbose);

    ifFeasible &= armourOptPtr1_->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, armourOptPtr1_->final_constr_violation);

    // second robot
    armourOptPtr2_->summarize_constraints(armourOptPtr2_->numCons, g + armourOptPtr1_->numCons, verbose);

    ifFeasible &= armourOptPtr2_->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, armourOptPtr2_->final_constr_violation);

    // arm-arm collision constraints
    const Index offset = armourOptPtr1_->numCons + armourOptPtr2_->numCons;
    const size_t num_arm_arm_collision = robotInfoPtr1_->tc_begin_and_end.size() * robotInfoPtr2_->tc_begin_and_end.size();

    for(Index i = 0; i < armourOptPtr1_->num_time_steps; i++){
        Index j = 0;
        for (Index j1 = 0; j1 < robotInfoPtr1_->tc_begin_and_end.size(); j1++) {
            for (Index j2 = 0; j2 < robotInfoPtr2_->tc_begin_and_end.size(); j2++) {
                if (g[i * num_arm_arm_collision + j + offset] < -constr_viol_tol) {
                    if (verbose) {
                        std::cout << "DualArmourOptimizer.cpp: "
                                  << " TC " << j1 << " on arm 1"
                                  << " and TC " << j2 << " on arm 2"
                                  << " collides at " << i << "th time step"
                                  << " with distance " << g[i * num_arm_arm_collision + j + offset] << "!\n";
                    }
                    final_constr_violation = std::max(final_constr_violation, -g[i * num_arm_arm_collision + j + offset]);
                }
                j++;
            }
        }
    }
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR