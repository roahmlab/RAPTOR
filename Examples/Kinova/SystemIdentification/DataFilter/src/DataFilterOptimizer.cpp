#include "DataFilterOptimizer.h"

namespace IDTO {
namespace Kinova {

bool DataFilterOptimizer::set_parameters(
        const VecX& x0_input,
        const VecX& tspan_input,
        const MatX& q_input,
        const MatX& q_d_input,
        const int degree_input,
        const int base_frequency_input) {
    enable_hessian = true;

    x0 = x0_input;

    if (q_input.rows() != q_d_input.rows() ||
        q_input.cols() != q_d_input.cols()) {
        THROW_EXCEPTION(IpoptException, "*** q_input and q_d_input should have the same size!");
    }

    if (q_input.cols() != tspan_input.size()) {
        THROW_EXCEPTION(IpoptException, "*** q_input and tspan_input should have the same size!");
    }

    q_data = q_input;
    q_d_data = q_d_input;

    // read joint limits from KinovaConstants.h
    VecX JOINT_LIMITS_LOWER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_LOWER, NUM_JOINTS);

    VecX JOINT_LIMITS_UPPER_VEC = Utils::initializeEigenVectorFromArray(JOINT_LIMITS_UPPER, NUM_JOINTS);

    // fixed frequency fourier curves with fixed initial position and velocity
    trajPtr_ = std::make_shared<FixedFrequencyFourierCurves>(tspan_input,
                                                             q_data.rows(), 
                                                             degree_input,
                                                             base_frequency_input,
                                                             q_data.col(0),
                                                             q_d_data.col(0));
    // trajPtr_ = std::make_shared<FourierCurves>(tspan_input,
    //                                            q_data.rows(), 
    //                                            degree_input,
    //                                            q_data.col(0),
    //                                            q_d_data.col(0));

    // // Joint limits
    // constraintsPtrVec_.push_back(std::make_unique<JointLimits>(trajPtr_, 
    //                                                            JOINT_LIMITS_LOWER_VEC, 
    //                                                            JOINT_LIMITS_UPPER_VEC));
    // constraintsNameVec_.push_back("joint limits");

    return true;
}

bool DataFilterOptimizer::get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
) {
    // number of decision variables
    numVars = trajPtr_->varLength;
    n = numVars;

    // number of inequality constraint
    numCons = 0;
    for ( Index i = 0; i < constraintsPtrVec_.size(); i++ ) {
        numCons += constraintsPtrVec_[i]->m;
    }
    m = numCons;

    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool DataFilterOptimizer::get_bounds_info(
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
        x_l[i] = -100;
    }

    // upper bounds  
    for( Index i = 0; i < n; i++ ) {
        x_u[i] = 100;
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

bool DataFilterOptimizer::eval_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number&       obj_value
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, false);

    // position estimation error
    Number position_error = 0;
    for ( Index i = 0; i < q_data.cols(); i++ ) {
        VecX qdiff = trajPtr_->q(i) - q_data.col(i);
        position_error += sqrt(qdiff.dot(qdiff));
    }

    // velocity estimation error
    Number velocity_error = 0;
    for ( Index i = 0; i < q_d_data.cols(); i++ ) {
        VecX q_d_diff = trajPtr_->q_d(i) - q_d_data.col(i);
        velocity_error += sqrt(q_d_diff.dot(q_d_diff));
    }

    obj_value = position_weight * position_error + 
                velocity_weight * velocity_error;

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool DataFilterOptimizer::eval_grad_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number*       grad_f
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, true);

    // gradient of position estimation error
    VecX grad_position_error = VecX::Zero(n);
    for ( Index i = 0; i < q_data.cols(); i++ ) {
        VecX q_diff = trajPtr_->q(i) - q_data.col(i);
        Number error = sqrt(q_diff.dot(q_diff));
        if (error > SQUARE_ROOT_THRESHOLD) {
            grad_position_error += q_diff.transpose() * trajPtr_->pq_pz(i) / error;
        }
        // else {
        //     grad_position_error += VecX::Zero(n);
        // }
    }

    // gradient of velocity estimation error
    VecX grad_velocity_error = VecX::Zero(n);
    for ( Index i = 0; i < q_d_data.cols(); i++ ) {
        VecX q_d_diff = trajPtr_->q_d(i) - q_d_data.col(i);
        Number error = sqrt(q_d_diff.dot(q_d_diff));
        if (error > SQUARE_ROOT_THRESHOLD) {
            grad_velocity_error += q_d_diff.transpose() * trajPtr_->pq_d_pz(i) / error;
        }
        // else {
        //     grad_velocity_error += VecX::Zero(n);
        // }
    }

    for ( Index i = 0; i < n; i++ ) {
        grad_f[i] = position_weight * grad_position_error(i) + 
                    velocity_weight * grad_velocity_error(i);
    }

    return true;
}

bool DataFilterOptimizer::eval_hess_f(
    Index         n,
    const Number* x,
    bool          new_x,
    MatX&         hess_f
) {
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    trajPtr_->compute(z, true, true);

    // hessian of position estimation error
    MatX hess_position_error = MatX::Zero(n, n);
    for ( Index i = 0; i < q_data.cols(); i++ ) {
        VecX q_diff = trajPtr_->q(i) - q_data.col(i);
        Number error = sqrt(q_diff.dot(q_diff));
        
        if (error > SQUARE_ROOT_THRESHOLD) {
            hess_position_error += trajPtr_->pq_pz(i).transpose() * trajPtr_->pq_pz(i) / error;

            for (Index j = 0; j < q_diff.size(); j++) {
                hess_position_error += q_diff(j) * trajPtr_->pq_pz_pz(j, i) / error;
            }

            MatX pdiffSquare_pz = q_diff.transpose() * trajPtr_->pq_pz(i);
            hess_position_error -= pdiffSquare_pz.transpose() * pdiffSquare_pz / std::pow(error, 3.0);
        }
        // else {
        //     hess_position_error += MatX::Zero(n, n);
        // }
    }

    // hessian of velocity estimation error
    MatX hess_velocity_error = MatX::Zero(n, n);
    for ( Index i = 0; i < q_d_data.cols(); i++ ) {
        VecX q_d_diff = trajPtr_->q_d(i) - q_d_data.col(i);
        Number error = sqrt(q_d_diff.dot(q_d_diff));
        
        if (error > SQUARE_ROOT_THRESHOLD) {
            hess_velocity_error += trajPtr_->pq_d_pz(i).transpose() * trajPtr_->pq_d_pz(i) / error;

            for (Index j = 0; j < q_d_diff.size(); j++) {
                hess_velocity_error += q_d_diff(j) * trajPtr_->pq_d_pz_pz(j, i) / error;
            }

            MatX pdiffSquare_pz = q_d_diff.transpose() * trajPtr_->pq_d_pz(i);
            hess_velocity_error -= pdiffSquare_pz.transpose() * pdiffSquare_pz / std::pow(error, 3.0);
        }
        // else {
        //     hess_position_error += MatX::Zero(n, n);
        // }
    }

    hess_f = position_weight * hess_position_error + 
             velocity_weight * hess_velocity_error;

    return true;
}

}; // namespace Kinova
}; // namespace IDTO