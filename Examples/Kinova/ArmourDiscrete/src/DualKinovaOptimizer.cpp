#include "DualKinovaOptimizer.h"

namespace RAPTOR {
namespace Kinova {

// // constructor
// DualKinovaOptimizer::DualKinovaOptimizer()
// {
// }


// // destructor
// DualKinovaOptimizer::~DualKinovaOptimizer()
// {
// }

// [TNLP_set_parameters]
bool DualKinovaOptimizer::set_parameters(
    const VecX& x0_input,
    const double T_input,
    const int N_input,
    const Model& model1_input, 
    const Model& model2_input, 
    const VecX& q0_input,
    const VecX& qT_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input,
    const VecX& joint_limits_buffer_input,
    const VecX& velocity_limits_buffer_input,
    const VecX& torque_limits_buffer_input,
    const bool include_gripper_or_not,
    const double collision_buffer_input,
    const double arm_arm_collision_buffer_input
 ) 
{
    enable_hessian = false;
    x0 = x0_input;
    arm_arm_collision_buffer = arm_arm_collision_buffer_input;

    kinovaOptPtr1_ = std::make_shared<KinovaLongerHorizonOptimizer>();
    kinovaOptPtr2_ = std::make_shared<KinovaLongerHorizonOptimizer>();

    kinovaOptPtr1_->set_parameters(x0_input.head(model1_input.nv), 
                                   T_input, 
                                   N_input, 
                                   degree_input, 
                                   model1_input, 
                                   q0_input.head(model1_input.nv), 
                                   qT_input.head(model1_input.nv),
                                   boxCenters_input, 
                                   boxOrientation_input, 
                                   boxSize_input, 
                                   joint_limits_buffer_input, 
                                   velocity_limits_buffer_input, 
                                   torque_limits_buffer_input, 
                                   include_gripper_or_not, 
                                   collision_buffer_input);

    kinovaOptPtr2_->set_parameters(x0_input.tail(model2_input.nv),
                                   T_input,
                                   N_input,
                                   degree_input,
                                   model2_input,
                                   q0_input.tail(model2_input.nv), 
                                   qT_input.tail(model2_input.nv),
                                   boxCenters_input,
                                   boxOrientation_input,
                                   boxSize_input,
                                   joint_limits_buffer_input,
                                   velocity_limits_buffer_input,
                                   torque_limits_buffer_input,
                                   include_gripper_or_not,
                                   collision_buffer_input);
                                                                                                                                                                                                                                                                                                                                                                        
    assert(x0.size() == kinovaOptPtr1_->trajPtr_->varLength + kinovaOptPtr2_->trajPtr_->varLength);

    return true;
}
// [TNLP_set_parameters]

// [TNLP_get_nlp_info]
// returns some info about the nlp
bool DualKinovaOptimizer::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // number of constraints
    numCons = 0;

    Index m1 = 0;
    kinovaOptPtr1_->get_nlp_info(n, m1, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m1;

    Index m2 = 0;
    kinovaOptPtr2_->get_nlp_info(n, m2, nnz_jac_g, nnz_h_lag, index_style);
    numCons += m2;

    // arm-arm collision constraints
    // for this example, we only check the last two tapered capsules of two arms
    const Index m_arm_arm = kinovaOptPtr1_->trajPtr_->N * 2 * 2;
    numCons += m_arm_arm;
    
    m = numCons;

    // number of decision variables
    numVars = kinovaOptPtr1_->numVars + kinovaOptPtr2_->numVars;
    n = numVars;

    // the nonzero structure of the Jacobian is the same for both arms, so they are stored separately
    nnz_jac_g = 
        kinovaOptPtr1_->numCons * kinovaOptPtr1_->trajPtr_->varLength + kinovaOptPtr2_->numCons * kinovaOptPtr2_->trajPtr_->varLength + 
        m_arm_arm * numVars;

    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool DualKinovaOptimizer::get_bounds_info(
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
    if(n != numVars){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_bounds_info!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in get_bounds_info!");
    }

    kinovaOptPtr1_->get_bounds_info(kinovaOptPtr1_->numVars, 
                                    x_l, 
                                    x_u, 
                                    kinovaOptPtr1_->numCons, 
                                    g_l, 
                                    g_u);
    kinovaOptPtr2_->get_bounds_info(kinovaOptPtr2_->numVars, 
                                    x_l + kinovaOptPtr1_->numVars, 
                                    x_u + kinovaOptPtr1_->numVars, 
                                    kinovaOptPtr2_->numCons, 
                                    g_l + kinovaOptPtr1_->numCons, 
                                    g_u + kinovaOptPtr1_->numCons);

    const Index offset = kinovaOptPtr1_->numCons + kinovaOptPtr1_->numCons;
    for (Index i = 0; i < kinovaOptPtr1_->trajPtr_->N * 2 * 2; i++) {
        g_l[offset + i] = arm_arm_collision_buffer;
        g_u[offset + i] = 1e19;
    }

    g_lb_copy = Utils::initializeEigenVectorFromArray(g_l, m);
    g_ub_copy = Utils::initializeEigenVectorFromArray(g_u, m);

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_eval_f]
// returns the value of the objective function
bool DualKinovaOptimizer::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    Number obj_value1 = 0;
    kinovaOptPtr1_->eval_f(kinovaOptPtr1_->numVars, x, new_x, obj_value1);

    Number obj_value2 = 0;
    kinovaOptPtr2_->eval_f(kinovaOptPtr2_->numVars, x + kinovaOptPtr1_->numVars, new_x, obj_value2);

    obj_value = obj_value1 + obj_value2;

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool DualKinovaOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    Number grad_f1[kinovaOptPtr1_->numVars];
    kinovaOptPtr1_->eval_grad_f(kinovaOptPtr1_->numVars, x, new_x, grad_f1);

    for(Index i = 0; i < kinovaOptPtr1_->numVars; i++){
        grad_f[i] = grad_f1[i];
    }

    Number grad_f2[kinovaOptPtr2_->numVars];
    kinovaOptPtr2_->eval_grad_f(kinovaOptPtr2_->numVars, x + kinovaOptPtr1_->numVars, new_x, grad_f2);

    for(Index i = 0; i < kinovaOptPtr2_->numVars; i++){
        grad_f[i + kinovaOptPtr1_->numVars] = grad_f2[i];
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool DualKinovaOptimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if(n != numVars){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_g!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "Error wrong value of m in eval_g!");
    }

    try {
        kinovaOptPtr1_->eval_g(kinovaOptPtr1_->numVars, x, new_x, kinovaOptPtr1_->numCons, g);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the first robot!");
    }

    try {
        kinovaOptPtr2_->eval_g(kinovaOptPtr2_->numVars, x + kinovaOptPtr1_->numVars, new_x, kinovaOptPtr2_->numCons, g + kinovaOptPtr1_->numCons);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g of the second robot!");
    }

    // arm-arm collision constraints
    const Index offset = kinovaOptPtr1_->numCons + kinovaOptPtr2_->numCons;
    const KinovaCustomizedConstraints* kccPtr1_ = dynamic_cast<
        const KinovaCustomizedConstraints*>(kinovaOptPtr1_->constraintsPtrVec_.back().get());
    const KinovaCustomizedConstraints* kccPtr2_ = dynamic_cast<
        const KinovaCustomizedConstraints*>(kinovaOptPtr2_->constraintsPtrVec_.back().get());
    const auto& sphere_radius1 = kccPtr1_->sphere_radius;
    const auto& sphere_radius2 = kccPtr2_->sphere_radius;
    const auto& tapered_capsules1 = kccPtr1_->tapered_capsules;
    const auto& tapered_capsules2 = kccPtr2_->tapered_capsules;

    try {
        for(Index i = 0; i < kinovaOptPtr1_->trajPtr_->N; i++){
            Index j = 0;
            for (Index j1 = 1; j1 < 3; j1++) {
                for (Index j2 = 1; j2 < 3; j2++) {
                    const size_t tc1_begin_index = tapered_capsules1[j1].first;
                    const size_t tc1_end_index   = tapered_capsules1[j1].second;
                    const size_t tc2_begin_index = tapered_capsules2[j2].first;
                    const size_t tc2_end_index   = tapered_capsules2[j2].second;

                    const Vec3& tc1_sphere_1 = kccPtr1_->sphere_centers_copy(tc1_begin_index, i);
                    const Vec3& tc1_sphere_2 = kccPtr1_->sphere_centers_copy(tc1_end_index, i);
                    const Vec3& tc2_sphere_1 = kccPtr2_->sphere_centers_copy(tc2_begin_index, i);
                    const Vec3& tc2_sphere_2 = kccPtr2_->sphere_centers_copy(tc2_end_index, i);

                    const double tc1_sphere_1_radius = sphere_radius1[tc1_begin_index];
                    const double tc1_sphere_2_radius = sphere_radius1[tc1_end_index];
                    const double tc2_sphere_1_radius = sphere_radius2[tc2_begin_index];
                    const double tc2_sphere_2_radius = sphere_radius2[tc2_end_index];

                    g[i * 4 + j + offset] = tcc.computeDistance(
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
bool DualKinovaOptimizer::eval_jac_g(
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
    if(n != numVars){
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
        for(Index i = 0; i < kinovaOptPtr1_->numCons; i++){
            for(Index j = 0; j < kinovaOptPtr1_->numVars; j++){
                iRow[idx] = i;
                jCol[idx] = j;
                idx++;
            }
        }

        // second robot
        for (Index i = 0; i < kinovaOptPtr2_->numCons; i++) {
            for (Index j = 0; j < kinovaOptPtr2_->numVars; j++) {
                iRow[idx] = i + kinovaOptPtr1_->numCons;
                jCol[idx] = j + kinovaOptPtr1_->numVars;
                idx++;
            }
        }

        // update iRow and jCol for arm-arm collision constraints
        Index offset = kinovaOptPtr1_->numCons + kinovaOptPtr2_->numCons;
        for (Index i = 0; i < kinovaOptPtr1_->trajPtr_->N; i++) {
            Index j = 0;
            for (Index j1 = 1; j1 < 3; j1++) {
                for (Index j2 = 1; j2 < 3; j2++) {
                    for (Index k = 0; k < numVars; k++) {
                        iRow[idx] = i * 4 + j + offset;
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
            kinovaOptPtr1_->eval_jac_g(kinovaOptPtr1_->numVars, x, new_x, kinovaOptPtr1_->numCons, nele_jac, iRow, jCol, values);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the first robot!");
        }

        try {
            kinovaOptPtr2_->eval_jac_g(kinovaOptPtr2_->numVars, x + kinovaOptPtr1_->numVars, new_x, kinovaOptPtr2_->numCons, nele_jac, iRow, jCol, values + kinovaOptPtr1_->numCons * kinovaOptPtr1_->numVars);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g of the second robot!");
        }

        // arm-arm collision constraints
        const Index offset = kinovaOptPtr1_->numCons * kinovaOptPtr1_->numVars + kinovaOptPtr2_->numCons * kinovaOptPtr2_->numVars;
        const KinovaCustomizedConstraints* kccPtr1_ = dynamic_cast<
            const KinovaCustomizedConstraints*>(kinovaOptPtr1_->constraintsPtrVec_.back().get());
        const KinovaCustomizedConstraints* kccPtr2_ = dynamic_cast<
            const KinovaCustomizedConstraints*>(kinovaOptPtr2_->constraintsPtrVec_.back().get());
        const auto& sphere_radius1 = kccPtr1_->sphere_radius;
        const auto& sphere_radius2 = kccPtr2_->sphere_radius;
        const auto& tapered_capsules1 = kccPtr1_->tapered_capsules;
        const auto& tapered_capsules2 = kccPtr2_->tapered_capsules;

        try {
            Eigen::Matrix<double, 3, numVars_dummy> dk_tc1_sphere_1;
            Eigen::Matrix<double, 3, numVars_dummy> dk_tc1_sphere_2;
            Eigen::Matrix<double, 3, numVars_dummy> dk_tc2_sphere_1;
            Eigen::Matrix<double, 3, numVars_dummy> dk_tc2_sphere_2;
            dk_tc1_sphere_1.setZero();
            dk_tc1_sphere_2.setZero();
            dk_tc2_sphere_1.setZero();
            dk_tc2_sphere_2.setZero();

            for(Index i = 0; i < kinovaOptPtr1_->trajPtr_->N; i++){
                Index j = 0;
                for (Index j1 = 1; j1 < 3; j1++) {
                    for (Index j2 = 1; j2 < 3; j2++) {
                        const size_t tc1_begin_index = tapered_capsules1[j1].first;
                        const size_t tc1_end_index   = tapered_capsules1[j1].second;
                        const size_t tc2_begin_index = tapered_capsules2[j2].first;
                        const size_t tc2_end_index   = tapered_capsules2[j2].second;

                        const Vec3& tc1_sphere_1 = kccPtr1_->sphere_centers_copy(tc1_begin_index, i);
                        const Vec3& tc1_sphere_2 = kccPtr1_->sphere_centers_copy(tc1_end_index, i);
                        const Vec3& tc2_sphere_1 = kccPtr2_->sphere_centers_copy(tc2_begin_index, i);
                        const Vec3& tc2_sphere_2 = kccPtr2_->sphere_centers_copy(tc2_end_index, i);

                        const double tc1_sphere_1_radius = sphere_radius1[tc1_begin_index];
                        const double tc1_sphere_2_radius = sphere_radius1[tc1_end_index];
                        const double tc2_sphere_1_radius = sphere_radius2[tc2_begin_index];
                        const double tc2_sphere_2_radius = sphere_radius2[tc2_end_index];

                        dk_tc1_sphere_1.leftCols(kinovaOptPtr1_->numVars) = kccPtr1_->sphere_centers_gradient_copy(tc1_begin_index, i);
                        dk_tc1_sphere_2.leftCols(kinovaOptPtr1_->numVars) = kccPtr1_->sphere_centers_gradient_copy(tc1_end_index, i);
                        dk_tc2_sphere_1.rightCols(kinovaOptPtr2_->numVars) = kccPtr2_->sphere_centers_gradient_copy(tc2_begin_index, i);
                        dk_tc2_sphere_2.rightCols(kinovaOptPtr2_->numVars) = kccPtr2_->sphere_centers_gradient_copy(tc2_end_index, i);

                        Eigen::Vector<double, numVars_dummy> dk_distances;
                        const double distance = tcc.computeDistance(
                            tc1_sphere_1, tc1_sphere_2, 
                            tc2_sphere_1, tc2_sphere_2,
                            dk_tc1_sphere_1, dk_tc1_sphere_2, 
                            dk_tc2_sphere_1, dk_tc2_sphere_2,
                            tc1_sphere_1_radius, tc1_sphere_2_radius, 
                            tc2_sphere_1_radius, tc2_sphere_2_radius,
                            dk_distances);

                        std::memcpy(values + (i * 4 + j) * numVars + offset, dk_distances.data(), numVars * sizeof(Number));

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

void DualKinovaOptimizer::summarize_constraints(
    Index                      m,
    const Number*              g,
    const bool                 verbose
) 
{
    ifFeasible = true;
    final_constr_violation = 0;

    // first robot
    kinovaOptPtr1_->summarize_constraints(kinovaOptPtr1_->numCons, g, verbose);

    ifFeasible &= kinovaOptPtr1_->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, kinovaOptPtr1_->final_constr_violation);

    // second robot
    kinovaOptPtr2_->summarize_constraints(kinovaOptPtr2_->numCons, g + kinovaOptPtr1_->numCons, verbose);

    ifFeasible &= kinovaOptPtr2_->ifFeasible;
    final_constr_violation = std::max(final_constr_violation, kinovaOptPtr2_->final_constr_violation);

    // arm-arm collision constraints
    const Index offset = kinovaOptPtr1_->numCons + kinovaOptPtr2_->numCons;
    for(Index i = 0; i < kinovaOptPtr1_->trajPtr_->N; i++){
        Index j = 0;
        for (Index j1 = 1; j1 < 3; j1++) {
            for (Index j2 = 1; j2 < 3; j2++) {
                if (g[i * 4 + j + offset] < -constr_viol_tol) {
                    if (verbose) {
                        std::cout << "DualKinovaOptimizer.cpp: "
                                  << " TC " << j1 << " on arm 1"
                                  << " and TC " << j2 << " on arm 2"
                                  << " collides at " << i << "th time step"
                                  << " with distance " << g[i * 4 + j + offset] << "!\n";
                    }
                    final_constr_violation = std::max(final_constr_violation, -g[i * 4 + j + offset]);
                }
                j++;
            }
        }
    }
}

}; // namespace Kinova
}; // namespace RAPTOR