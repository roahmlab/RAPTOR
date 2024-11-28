#include "ArmourOptimizer.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

bool ArmourOptimizer::set_parameters(
    const VecX& q_des_input,
    Number t_plan_input,
    const std::shared_ptr<RobotInfo>& robotInfoPtr_input,
    const std::shared_ptr<BezierCurveInterval>& trajPtr_input,
    const std::shared_ptr<PZDynamics>& dynPtr_input,
    const std::vector<Vec3>& boxCenters_input,
    const std::vector<Vec3>& boxOrientation_input,
    const std::vector<Vec3>& boxSize_input
 ) 
 {
    enable_hessian = false;

    q_des = q_des_input;
    t_plan = t_plan_input;

    robotInfoPtr_ = robotInfoPtr_input;
    num_spheres = robotInfoPtr_->num_spheres;
    num_fixed_joints = robotInfoPtr_->num_joints - robotInfoPtr_->num_motors;

    if (num_fixed_joints < 0 || num_fixed_joints > 1) {
        THROW_EXCEPTION(IpoptException, "*** Error: num_fixed_joints is not 0 or 1!");
    }

    trajPtr_ = trajPtr_input;
    num_time_steps = trajPtr_->num_time_steps;

    dynPtr_ = dynPtr_input;

    bcaPtrs.resize(num_time_steps);
    for (size_t i = 0; i < num_time_steps; i++) {
        bcaPtrs[i] = std::make_shared<BoxCollisionAvoidance>(
            boxCenters_input, boxOrientation_input, boxSize_input);
        bcaPtrs[i]->onlyComputeDerivativesForMinimumDistance = true;
    }

    return true;
}

bool ArmourOptimizer::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
    // The problem described NUM_FACTORS variables, x[NUM_FACTORS] through x[NUM_FACTORS] for each joint
    numVars = NUM_FACTORS;
    n = NUM_FACTORS;

    // number of inequality constraint
    numCons = NUM_FACTORS * num_time_steps + // torque limits
              NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps + // contact constraints
              num_time_steps * num_spheres + // obstacle avoidance constraints
              NUM_FACTORS * 4; // joint position, velocity limits
    m = numCons;

    nnz_jac_g = m * n;
    nnz_h_lag = n * n;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool ArmourOptimizer::get_bounds_info(
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
    if(n != NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in get_bounds_info!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "*** Error wrong value of m in get_bounds_info!");
    }

    // lower bounds
    for( Index i = 0; i < n; i++ ) {
        x_l[i] = -1.0;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 1.0;
    }

    Index offset = 0;
    
    // torque limits
    for( Index i = 0; i < num_time_steps; i++ ) {
        for( Index j = 0; j < NUM_FACTORS; j++ ) {
            g_l[i * NUM_FACTORS + j] = TORQUE_LIMITS_LOWER[j] + dynPtr_->torque_radii(j, i);
            g_u[i * NUM_FACTORS + j] = TORQUE_LIMITS_UPPER[j] - dynPtr_->torque_radii(j, i);
        }
    }    
    offset += NUM_FACTORS * num_time_steps;

    // contact constraints
    for( Index i = offset; i < offset + NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps; i++ ) {
        g_l[i] = 0.0;
        g_u[i] = 1e19;
    }
    offset += NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps;

    // collision avoidance constraints
    for( Index i = offset; i < offset + num_time_steps * num_spheres; i++ ) {
        g_l[i] = 0.0;
        g_u[i] = 1e19;
    }
    offset += num_time_steps * num_spheres;

    // state limit constraints
    //     minimum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = Utils::deg2rad(JOINT_LIMITS_LOWER[i - offset]) + robotInfoPtr_->ultimate_bound_info.qe;
        g_u[i] = Utils::deg2rad(JOINT_LIMITS_UPPER[i - offset]) - robotInfoPtr_->ultimate_bound_info.qe;
    }
    offset += NUM_FACTORS;

    //     maximum joint position
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = Utils::deg2rad(JOINT_LIMITS_LOWER[i - offset]) + robotInfoPtr_->ultimate_bound_info.qe;
        g_u[i] = Utils::deg2rad(JOINT_LIMITS_UPPER[i - offset]) - robotInfoPtr_->ultimate_bound_info.qe;
    }
    offset += NUM_FACTORS;

    //     minimum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = Utils::deg2rad(VELOCITY_LIMITS_LOWER[i - offset]) + robotInfoPtr_->ultimate_bound_info.qde;
        g_u[i] = Utils::deg2rad(VELOCITY_LIMITS_UPPER[i - offset]) - robotInfoPtr_->ultimate_bound_info.qde;
    }
    offset += NUM_FACTORS;

    //     maximum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        g_l[i] = Utils::deg2rad(VELOCITY_LIMITS_LOWER[i - offset]) + robotInfoPtr_->ultimate_bound_info.qde;
        g_u[i] = Utils::deg2rad(VELOCITY_LIMITS_UPPER[i - offset]) - robotInfoPtr_->ultimate_bound_info.qde;
    }

    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool ArmourOptimizer::get_starting_point(
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

    if(n != NUM_FACTORS){
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
bool ArmourOptimizer::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
    if(n != NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_f!");
    }

    // obj_value = sum((q_plan - q_des).^2);
    VecX q_plan(n);
    for(Index i = 0; i < n; i++){
        const Number k_actual = trajPtr_->k_center[i] + trajPtr_->k_range[i] * x[i];
        q_plan[i] = q_des_func(
            trajPtr_->q0[i], trajPtr_->Tqd0[i], trajPtr_->TTqdd0[i], 
            k_actual, t_plan);
    }

    // kinova has 4 infinite rotation joints
    obj_value = pow(Utils::wrapToPi(q_plan[0] - q_des[0]), 2) + // These are continuous joints
                pow(Utils::wrapToPi(q_plan[2] - q_des[2]), 2) + 
                pow(Utils::wrapToPi(q_plan[4] - q_des[4]), 2) + 
                pow(Utils::wrapToPi(q_plan[6] - q_des[6]), 2) + 
                pow(q_plan[1] - q_des[1], 2) +           // These are not continuous joints
                pow(q_plan[3] - q_des[3], 2) + 
                pow(q_plan[5] - q_des[5], 2);

    update_minimal_cost_solution(n, Utils::initializeEigenVectorFromArray(x, n), new_x, obj_value);

    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool ArmourOptimizer::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    if(n != NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_grad_f!");
    }

    for (Index i = 0; i < n; i++) {
        const Number k_actual = trajPtr_->k_center[i] + trajPtr_->k_range[i] * x[i];
        Number q_plan = q_des_func(
            trajPtr_->q0[i], trajPtr_->Tqd0[i], trajPtr_->TTqdd0[i], 
            k_actual, t_plan);
        Number dk_q_plan = pow(t_plan,3) * (6 * pow(t_plan,2) - 15 * t_plan + 10) * trajPtr_->k_range[i];

        // kinova has 4 infinite rotation joints
        if (i % 2 == 0) {
            grad_f[i] = (2 * Utils::wrapToPi(q_plan - q_des[i]) * dk_q_plan);
        }
        else {
            grad_f[i] = (2 * (q_plan - q_des[i]) * dk_q_plan);
        }
    }

    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool ArmourOptimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if(n != NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_g!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "Error wrong value of m in eval_g!");
    }

    try {
        Index i = 0, offset = 0;

        // torque limits
        try {
            #pragma omp parallel for shared(dynPtr_, x, g) private(i) schedule(dynamic)
            for (i = 0; i < num_time_steps; i++) {
                for (size_t j = 0; j < NUM_FACTORS; j++) {
                    const PZSparse& PZtorque = dynPtr_->data_sparses[i].tau(j);
                    const Interval res = PZtorque.slice(x);
                    g[i * NUM_FACTORS + j] = getCenter(res);
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_g!");
        }                                                                                                                                                                           offset += num_time_steps * NUM_FACTORS;
        offset += num_time_steps * NUM_FACTORS;
        
        // contact constraints (separation, friction cone, ZMP)
        if (num_fixed_joints > 0) {
            try {
                #pragma omp parallel for shared(dynPtr_, x, g) private(i) schedule(dynamic)
                for (i = 0; i < num_time_steps; i++) {
                    // (1) support force larger than 0
                    const Interval supportForceRange = 
                        dynPtr_->data_sparses[i].f[dynPtr_->model_sparses[i].nv].linear()(2)
                            .slice(x);
                    g[i * NUM_CONTACT_CONSTRAINTS + offset] = supportForceRange.lower();

                    // (2) friction cone constraints
                    for (size_t j = 0; j < FRICTION_CONE_LINEARIZED_SIZE; j++) {
                        const Interval frictionConstraintRange = dynPtr_->friction_PZs(j, i).slice(x);
                        g[i * NUM_CONTACT_CONSTRAINTS + 1 + j + offset] = frictionConstraintRange.lower();
                    }      

                    // (3) ZMP constraints
                    for (size_t j = 0; j < ZMP_LINEARIZED_SIZE; j++) {
                        const Interval zmpConstraintRange = dynPtr_->zmp_PZs(j, i).slice(x);
                        g[i * NUM_CONTACT_CONSTRAINTS + 1 + FRICTION_CONE_LINEARIZED_SIZE + j + offset] = zmpConstraintRange.lower();
                    }
                }
            }
            catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "Error in eval_g!");
            }
            offset += NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps;
        }

        // obstacle avoidance constraints
        try {
            #pragma omp parallel for shared(dynPtr_, bcaPtrs, x, g) private(i) schedule(dynamic)
            for (i = 0; i < num_time_steps; i++) {
                for (size_t j = 0; j < num_spheres; j++) {
                    const std::string sphere_name = "collision-" + std::to_string(j);
                    const pinocchio::FrameIndex frame_id = 
                        robotInfoPtr_->model.getFrameId(sphere_name);

                    const auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();

                    const Interval x_res = PZsphere(0).slice(x);
                    const Interval y_res = PZsphere(1).slice(x);
                    const Interval z_res = PZsphere(2).slice(x);

                    Vec3 sphere_center;
                    sphere_center << getCenter(x_res), 
                                     getCenter(y_res), 
                                     getCenter(z_res);

                    bcaPtrs[i]->computeDistance(sphere_center);
                    g[i * num_spheres + j + offset] = bcaPtrs[i]->minimumDistance - dynPtr_->sphere_radii(j, i);
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_g!");
        }
        offset += num_time_steps * num_spheres;

        trajPtr_->returnJointPositionExtremum(g + offset, x);
        offset += NUM_FACTORS * 2;

        trajPtr_->returnJointVelocityExtremum(g + offset, x);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        THROW_EXCEPTION(IpoptException, "Error in eval_g!");
    }

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool ArmourOptimizer::eval_jac_g(
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
    if(n != NUM_FACTORS){
        THROW_EXCEPTION(IpoptException, "Error wrong value of n in eval_jac_g!");
    }
    if(m != numCons){
        THROW_EXCEPTION(IpoptException, "Error wrong value of m in eval_jac_g!");
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
        Index i = 0, offset = 0;

        // torque limits
        try {
            #pragma omp parallel for shared(dynPtr_, x, values) private(i) schedule(dynamic)
            for(i = 0; i < num_time_steps; i++) {
                for (size_t j = 0; j < NUM_FACTORS; j++) {
                    const PZSparse& PZtorque = dynPtr_->data_sparses[i].tau(j);
                    PZtorque.slice(values + (i * NUM_FACTORS + j) * NUM_FACTORS, x);
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g!");
        }
        offset += num_time_steps * NUM_FACTORS * NUM_FACTORS;

        // contact constraints gradient
        if (num_fixed_joints > 0) {
            try {
                #pragma omp parallel for shared(dynPtr_, x, values) private(i) schedule(dynamic)
                for (i = 0; i < num_time_steps; i++) {
                    // separation
                    const auto& PZseparation = dynPtr_->data_sparses[i].f[dynPtr_->model_sparses[i].nv].linear()(2);
                    PZseparation.slice(values + (i * NUM_CONTACT_CONSTRAINTS * num_fixed_joints) * NUM_FACTORS + offset, x);

                    // friction cone
                    for (size_t j = 0; j < FRICTION_CONE_LINEARIZED_SIZE; j++) {
                        const auto& PZfriction = dynPtr_->friction_PZs(j, i);
                        PZfriction.slice(values + (i * NUM_CONTACT_CONSTRAINTS * num_fixed_joints + 1 + j) * NUM_FACTORS + offset, x);
                    }

                    // ZMP
                    for (size_t j = 0; j < ZMP_LINEARIZED_SIZE; j++) {
                        const auto& PZzmp = dynPtr_->zmp_PZs(j, i);
                        PZzmp.slice(values + (i * NUM_CONTACT_CONSTRAINTS * num_fixed_joints + 1 + FRICTION_CONE_LINEARIZED_SIZE + j) * NUM_FACTORS + offset, x);
                    }
                }
            }
            catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
                THROW_EXCEPTION(IpoptException, "Error in eval_jac_g!");
            }
            offset += NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps * NUM_FACTORS;
        }

        // obstacle avoidance constraints gradient
        try {
            #pragma omp parallel for shared(dynPtr_, x, values) private(i) schedule(dynamic)
            for(i = 0; i < num_time_steps; i++) {
                for (size_t j = 0; j < num_spheres; j++) {
                    const std::string sphere_name = "collision-" + std::to_string(j);
                    const pinocchio::FrameIndex frame_id = 
                        robotInfoPtr_->model.getFrameId(sphere_name);

                    const auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();

                    const Interval x_res = PZsphere(0).slice(x);
                    const Interval y_res = PZsphere(1).slice(x);
                    const Interval z_res = PZsphere(2).slice(x);

                    Vec3 sphere_center;
                    sphere_center << getCenter(x_res), 
                                     getCenter(y_res), 
                                     getCenter(z_res);

                    VecX dk_x_res(NUM_FACTORS), dk_y_res(NUM_FACTORS), dk_z_res(NUM_FACTORS);
                    PZsphere(0).slice(dk_x_res, x);
                    PZsphere(1).slice(dk_y_res, x);
                    PZsphere(2).slice(dk_z_res, x);

                    MatX dk_sphere_center(3, NUM_FACTORS);
                    dk_sphere_center.row(0) = dk_x_res;
                    dk_sphere_center.row(1) = dk_y_res;
                    dk_sphere_center.row(2) = dk_z_res;

                    bcaPtrs[i]->computeDistance(sphere_center, dk_sphere_center);

                    const VecX& dk_distances = bcaPtrs[i]->pdistances_pz.row(bcaPtrs[i]->minimumDistanceIndex);
                    std::memcpy(values + (i * num_spheres + j) * NUM_FACTORS + offset, dk_distances.data(), NUM_FACTORS * sizeof(Number));
                }
            }
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            THROW_EXCEPTION(IpoptException, "Error in eval_jac_g!");
        }
        offset += num_time_steps * num_spheres * NUM_FACTORS;

        // state limit constraints
        trajPtr_->returnJointPositionExtremumGradient(values + offset, x);
        offset += NUM_FACTORS * 2 * NUM_FACTORS;

        trajPtr_->returnJointVelocityExtremumGradient(values + offset, x);
    }

    return true;
}
// [TNLP_eval_jac_g]

void ArmourOptimizer::summarize_constraints(
    Index                      m,
    const Number*              g,
    const bool                 verbose
) 
{
    ifFeasible = true;
    final_constr_violation = 0;
    Index offset = 0;

    // control input constraints
    for( Index i = 0; i < num_time_steps; i++ ) {
        for( Index j = 0; j < NUM_FACTORS; j++ ) {
            const Number constr_violation = fmax(
                g[i * NUM_FACTORS + j] - (TORQUE_LIMITS_UPPER[j] - dynPtr_->torque_radii(j, i)), 
                (TORQUE_LIMITS_LOWER[j] + dynPtr_->torque_radii(j, i)) - g[i * NUM_FACTORS + j]);

            if (constr_violation > final_constr_violation) {
                final_constr_violation = constr_violation;
            }

            if (constr_violation > constr_viol_tol) {
                ifFeasible = false;
                
                if (verbose) {
                    std::cout << "ArmourOptimizer.cpp: Control torque of motor " << j + 1 << 
                                " at time interval " << i << " exceeds limit!\n";
                    std::cout << "    value: " << g[i * NUM_FACTORS + j] << "\n";
                    std::cout << "    range: [ " << TORQUE_LIMITS_LOWER[j] + dynPtr_->torque_radii(j, i) << ", "
                                                 << TORQUE_LIMITS_UPPER[j] - dynPtr_->torque_radii(j, i) << " ]\n";
                }
            }
        }
    }    
    offset += NUM_FACTORS * num_time_steps;

    // contact constraints validation
    if (num_fixed_joints > 0) {
        for( Index i = 0; i < num_time_steps; i++ ) {
            for( Index j = 0; j < NUM_CONTACT_CONSTRAINTS; j++ ) {
                const Number constr_violation = -g[i * NUM_CONTACT_CONSTRAINTS + j + offset];

                if (constr_violation > final_constr_violation) {
                    final_constr_violation = constr_violation;
                }

                if (constr_violation > constr_viol_tol) {
                    ifFeasible = false;
                    
                    if (verbose) {
                        if (j < 1) {
                            std::cout << "ArmourOptimizer.cpp: Separation constraint " << j << 
                                         " at time interval " << i << " is violated!\n";
                        }
                        else if (j < 1 + FRICTION_CONE_LINEARIZED_SIZE) {
                            std::cout << "ArmourOptimizer.cpp: Friction cone constraint " << j - 1 << 
                                         " at time interval " << i << " is violated!\n";
                        }
                        else {
                            std::cout << "ArmourOptimizer.cpp: ZMP constraint " << j - (1 + FRICTION_CONE_LINEARIZED_SIZE) << 
                                         " at time interval " << i << " is violated!\n";
                        }
                        std::cout << "    value: " << g[i * NUM_CONTACT_CONSTRAINTS + j + offset] << "\n";
                    }
                }
            } 
        }

        offset += NUM_CONTACT_CONSTRAINTS * num_fixed_joints * num_time_steps;
    }

    // obstacle avoidance constraints
    for( Index i = 0; i < num_time_steps; i++ ) {
        for( Index j = 0; j < num_spheres; j++ ) {
            const Number constr_violation = -g[i * num_spheres + j + offset];

            if (constr_violation > final_constr_violation) {
                final_constr_violation = constr_violation;
            }

            if (constr_violation > constr_viol_tol) {
                ifFeasible = false;
                
                if (verbose) {
                    std::cout << "ArmourOptimizer.cpp: Sphere " << j << 
                                " at time interval " << i << " collides with obstacles!\n";
                    std::cout << "    distance: " << g[i * num_spheres + j + offset] << "\n";
                }
            }
        }
    }
    offset += num_time_steps * num_spheres;

    // state limit constraints
    //     minimum joint position
    for( Index i = offset; i < offset + 2 * NUM_FACTORS; i++ ) {
        const Number constr_violation = fmax(
            g[i] - (Utils::deg2rad(JOINT_LIMITS_UPPER[(i - offset) % NUM_FACTORS]) - robotInfoPtr_->ultimate_bound_info.qe), 
            (Utils::deg2rad(JOINT_LIMITS_LOWER[(i - offset) % NUM_FACTORS]) + robotInfoPtr_->ultimate_bound_info.qe) - g[i]);

        if (constr_violation > final_constr_violation) {
            final_constr_violation = constr_violation;
        }

        if (constr_violation > constr_viol_tol) {
            ifFeasible = false;
                
            if (verbose) {
                std::cout << "ArmourOptimizer.cpp: joint " << i + 1 - offset << " exceeds position limit when it reaches minimum!\n";
                std::cout << "    value: " << g[i] << "\n";
                std::cout << "    range: [ " << Utils::deg2rad(JOINT_LIMITS_LOWER[i - offset]) + 
                                                robotInfoPtr_->ultimate_bound_info.qe << ", "
                                             << Utils::deg2rad(JOINT_LIMITS_UPPER[i - offset]) - 
                                                robotInfoPtr_->ultimate_bound_info.qe << " ]\n";
            }
        }
    }
    offset += 2 * NUM_FACTORS;

    //     minimum joint velocity
    for( Index i = offset; i < offset + NUM_FACTORS; i++ ) {
        const Number constr_violation = fmax(
            g[i] - (Utils::deg2rad(VELOCITY_LIMITS_UPPER[(i - offset) % NUM_FACTORS]) - robotInfoPtr_->ultimate_bound_info.qde), 
            (Utils::deg2rad(VELOCITY_LIMITS_LOWER[(i - offset) % NUM_FACTORS]) - robotInfoPtr_->ultimate_bound_info.qde) - g[i]);

        if (constr_violation > final_constr_violation) {
            final_constr_violation = constr_violation;
        }

        if (constr_violation > constr_viol_tol) {
            ifFeasible = false;
                
            if (verbose) {
                std::cout << "ArmourOptimizer.cpp: joint " << i + 1 - offset << " exceeds velocity limit when it reaches minimum!\n";
                std::cout << "    value: " << g[i] << "\n";
                std::cout << "    range: [ " << Utils::deg2rad(VELOCITY_LIMITS_LOWER[i - offset]) + 
                                                robotInfoPtr_->ultimate_bound_info.qde << ", "
                                             << Utils::deg2rad(VELOCITY_LIMITS_UPPER[i - offset]) - 
                                                robotInfoPtr_->ultimate_bound_info.qde << " ]\n";
            }
        }
    }
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR