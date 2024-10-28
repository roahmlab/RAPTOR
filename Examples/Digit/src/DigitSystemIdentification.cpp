#include "DigitSystemIdentification.h"

namespace RAPTOR {
namespace DigitWholeBodySysID {

// // constructor
// DigitSystemIdentification::DigitSystemIdentification()
// {
// }

// // destructors
// DigitSystemIdentification::~DigitSystemIdentification()
// {
// }

bool DigitSystemIdentification::set_parameters(
    const Model& model_input,
    const std::shared_ptr<MatX>& posPtr_input,
    const std::shared_ptr<MatX>& torquePtr_input
)
{ 
    enable_hessian = false;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    ddcPtr_ = std::make_shared<DigitWholeBodyDynamicsConstraints>(modelPtr_);

    Nact = NUM_INDEPENDENT_JOINTS;

    posPtr_ = posPtr_input;
    torquePtr_ = torquePtr_input;

    if (posPtr_->rows() != modelPtr_->nv) {
        throw std::invalid_argument("Error: input data matrices have different number of rows");
    }

    if (torquePtr_->rows() != Nact) {
        throw std::invalid_argument("Error: torque data matrix has different number of rows");
    }

    if (posPtr_->cols() != torquePtr_->cols()) {
        throw std::invalid_argument("Error: input data matrices have different number of columns");
    }

    N = posPtr_->cols();

    // find links that are not trivial (non-zero mass and inertia)
    nontrivialLinkIds.clear();
    // for (int i = 0; i < modelPtr_->nv; i++) {
    //     const int pinocchio_joint_id = i + 1;

    //     if (modelPtr_->inertias[pinocchio_joint_id].mass() > 0.0) {
    //         nontrivialLinkIds.push_back(i);
    //     }
    // }
    nontrivialLinkIds.push_back(modelPtr_->getJointId("Rz") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("left_hip_roll") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("left_hip_yaw") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("left_hip_pitch") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_ach2") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("left_knee") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_tarsus") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_toe_A") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_A2") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_toe_B") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_B2") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_toe_pitch") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("left_toe_roll") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("Rz") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("right_hip_roll") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("right_hip_yaw") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("right_hip_pitch") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_ach2") - 1);
    nontrivialLinkIds.push_back(modelPtr_->getJointId("right_knee") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_tarsus") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_toe_A") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_A2") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_toe_B") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_B2") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_toe_pitch") - 1);
    // nontrivialLinkIds.push_back(modelPtr_->getJointId("right_toe_roll") - 1);

    // directly set up initial condition here
    int n = 4 * nontrivialLinkIds.size() + 
            (NUM_DEPENDENT_JOINTS + 6) * N;

    x0 = VecX::Zero(n);
    for (int i = 0; i < nontrivialLinkIds.size(); i++) {
        const int pinocchio_joint_id = nontrivialLinkIds[i] + 1;
        const VecX& dynParams = modelPtr_->inertias[pinocchio_joint_id].toDynamicParameters();
        x0.segment(4 * i, 4) = dynParams.head(4);
    }

    const int CONSTRAINT_DIM = NUM_DEPENDENT_JOINTS + 6;

    A = MatX::Zero(NUM_INDEPENDENT_JOINTS * N, 4 * nontrivialLinkIds.size() + CONSTRAINT_DIM * N);
    b = VecX::Zero(NUM_INDEPENDENT_JOINTS * N);

    MatX A2 = MatX::Zero(NUM_DEPENDENT_JOINTS * N, 4 * nontrivialLinkIds.size() + CONSTRAINT_DIM * N);
    VecX b2 = VecX::Zero(NUM_DEPENDENT_JOINTS * N);

    for (int i = 0; i < N; i++) {
        const VecX& q = posPtr_->col(i);
        const VecX& tau = torquePtr_->col(i);

        pinocchio::computeJointTorqueRegressor(
            *modelPtr_, *dataPtr_, 
            q, VecX::Zero(modelPtr_->nv), VecX::Zero(modelPtr_->nv));

        for (int j = 0; j < nontrivialLinkIds.size(); j++) {
            MatX jointTorqueRegressorIndep(NUM_INDEPENDENT_JOINTS, 10);
            MatX jointTorqueRegressorDep(NUM_DEPENDENT_JOINTS, 10);

            ddcPtr_->get_independent_rows(
                jointTorqueRegressorIndep, 
                dataPtr_->jointTorqueRegressor.middleCols(10 * nontrivialLinkIds[j], 4));
            ddcPtr_->get_dependent_rows(
                jointTorqueRegressorDep, 
                dataPtr_->jointTorqueRegressor.middleCols(10 * nontrivialLinkIds[j], 4));

            A.block(i * NUM_INDEPENDENT_JOINTS, 
                    j * 4, 
                    NUM_INDEPENDENT_JOINTS, 
                    4) = jointTorqueRegressorIndep;
            A2.block(i * NUM_DEPENDENT_JOINTS, 
                     j * 4,
                     NUM_DEPENDENT_JOINTS, 
                     4) = jointTorqueRegressorDep;
        }

        tau_fixed = VecX::Zero(modelPtr_->nv);
        for (int j = 0; j < modelPtr_->nv; j++) {
            if (std::find(nontrivialLinkIds.begin(), nontrivialLinkIds.end(), j) != nontrivialLinkIds.end()) {
                continue;
            }

            const int pinocchio_joint_id = j + 1;
            tau_fixed += dataPtr_->jointTorqueRegressor.middleCols(10 * j, 10) * 
                modelPtr_->inertias[pinocchio_joint_id].toDynamicParameters();
        }

        ddcPtr_->get_J(q);

        MatX J_indep(CONSTRAINT_DIM, NUM_INDEPENDENT_JOINTS);
        MatX J_dep(CONSTRAINT_DIM, NUM_DEPENDENT_JOINTS);
        ddcPtr_->get_independent_columns(J_indep, ddcPtr_->J);
        ddcPtr_->get_dependent_columns(J_dep, ddcPtr_->J);

        A.block(i * NUM_INDEPENDENT_JOINTS, 
                4 * nontrivialLinkIds.size() + i * CONSTRAINT_DIM, 
                NUM_INDEPENDENT_JOINTS, 
                CONSTRAINT_DIM) = -J_indep.transpose();
        A2.block(i * NUM_DEPENDENT_JOINTS, 
                 4 * nontrivialLinkIds.size() + i * CONSTRAINT_DIM,
                 NUM_DEPENDENT_JOINTS, 
                 CONSTRAINT_DIM) = -J_dep.transpose();

        b.segment(i * NUM_INDEPENDENT_JOINTS, NUM_INDEPENDENT_JOINTS) = 
            tau - ddcPtr_->get_independent_vector(tau_fixed);
        b2.segment(i * NUM_DEPENDENT_JOINTS, NUM_DEPENDENT_JOINTS) =
            -ddcPtr_->get_dependent_vector(tau_fixed);

        Utils::writeEigenMatrixToFile(A, "A.txt");
        Utils::writeEigenMatrixToFile(A2, "A2.txt");

        // // verification
        // pinocchio::rnea(*modelPtr_, *dataPtr_, q, VecX::Zero(modelPtr_->nv), VecX::Zero(modelPtr_->nv));
        // VecX tau_verify(modelPtr_->nv);
        // ddcPtr_->fill_independent_vector(tau_verify, A.block(i * NUM_INDEPENDENT_JOINTS, 0, NUM_INDEPENDENT_JOINTS, 4 * nontrivialLinkIds.size()) * x0.head(4 * nontrivialLinkIds.size()));
        // ddcPtr_->fill_dependent_vector(tau_verify, A2.block(i * NUM_DEPENDENT_JOINTS, 0, NUM_DEPENDENT_JOINTS, 4 * nontrivialLinkIds.size()) * x0.head(4 * nontrivialLinkIds.size()));
        // tau_verify += tau_fixed;
        // std::cout << "Verification: " << (dataPtr_->tau - tau_verify).norm() << std::endl;
    }

    // // initialize LMI constraints for all links
    // constraintsPtrVec_.push_back(std::make_unique<LMIConstraints>(nontrivialLinkIds.size(), n));       
    // constraintsNameVec_.push_back("LMI constraints"); 

    // net force of the system to be 0
    const double totalMass = pinocchio::computeTotalMass(*modelPtr_);
    constraintsPtrVec_.push_back(std::make_unique<CustomizedConstraints>(totalMass, N, 4 * nontrivialLinkIds.size(), n));
    constraintsNameVec_.push_back("Net force of the system");

    // unconstrained part of the dynamics needs to be zero
    constraintsPtrVec_.push_back(std::make_unique<LinearConstraints>(A2, b2));
    constraintsNameVec_.push_back("Unconstrained part of the dynamics");

    return true;
}

bool DigitSystemIdentification::get_nlp_info(
    Index &n,
    Index &m,
    Index &nnz_jac_g,
    Index &nnz_h_lag,
    IndexStyleEnum &index_style
)
{
    // number of decision variables
    n = 4 * nontrivialLinkIds.size() + (NUM_DEPENDENT_JOINTS + 6) * N;
    numVars= n;

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

bool DigitSystemIdentification::get_bounds_info(
    Index n, 
    Number* x_l, 
    Number* x_u,
    Index m, 
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

    // use base class function to set bounds in g_l and g_u
    Optimizer::get_bounds_info(n, x_l, x_u, m, g_l, g_u);

    // set variable bounds (overwrite previous bounds in x_l and x_u)
    for (Index i = 0; i < 4 * nontrivialLinkIds.size(); i++) {
        if (i == 3) {
            x_l[i] = 0;
            x_u[i] = x0(i) + 0.5;
        }
        else {
            if (std::abs(x0(i)) > 0.1) {
                if (x0(i) > 0) {
                    x_l[i] = (1 - default_maximum_uncertainty) * x0(i);
                    x_u[i] = (1 + default_maximum_uncertainty) * x0(i);
                }
                else {
                    x_l[i] = (1 + default_maximum_uncertainty) * x0(i);
                    x_u[i] = (1 - default_maximum_uncertainty) * x0(i);
                }
            }
            else {
                x_l[i] = x0(i) - default_maximum_uncertainty;
                x_u[i] = x0(i) + default_maximum_uncertainty;
            }

            // mass >= 0
            if (i % 4 == 0) {
                x_l[i] = std::max(x_l[i], 0.0);
            }
        }
    }

        // ground reaction forces
    for (Index i = 10 * nontrivialLinkIds.size(); i < n; i++) {
        x_l[i] = -1e19;
        x_u[i] = 1e19;
    }

    return true;
}

bool DigitSystemIdentification::eval_f(
    Index n,
    const Number *x,
    bool new_x,
    Number &obj_value
)
{
    if(n != numVars){
       THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
    }

    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    VecX diff = A * z - b;

    obj_value = 0.5 * diff.squaredNorm() / N;

    update_minimal_cost_solution(n, z, new_x, obj_value);

    return true;
}

bool DigitSystemIdentification::eval_grad_f(
    Index n,
    const Number *x,
    bool new_x,
    Number *grad_f
)
{
    VecX z = Utils::initializeEigenVectorFromArray(x, n);
    VecX grad_f_vec = VecX::Zero(n);

    VecX diff = A * z - b;
    grad_f_vec = A.transpose() * diff;

    for (Index i = 0; i < n; i++) {
        grad_f[i] = grad_f_vec(i) / N;
    }

    return true;
}

bool DigitSystemIdentification::eval_hess_f(
    Index         n,
    const Number* x,
    bool          new_x,
    MatX&         hess_f
) {
    VecX z = Utils::initializeEigenVectorFromArray(x, n);

    hess_f = A.transpose() * A / N;

    return true;
}

void CustomizedConstraints::compute(const VecX& z, 
                                    bool compute_derivatives,
                                    bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    const VecX& lambdas = z.tail(N * (NUM_DEPENDENT_JOINTS + 6));

    for (int i = 0; i < N; i++) {
        const VecX& lambda = lambdas.segment(i * (NUM_DEPENDENT_JOINTS + 6), NUM_DEPENDENT_JOINTS + 6);
        const VecX& lambda_contact = lambda.tail(12);

        g.segment(3 * i, 3) = lambda_contact.head(3) + 
                              lambda_contact.segment(6, 3) + 
                              totalMass * Eigen::Vector3d(0, 0, -9.806);
        
        if (compute_derivatives) {
            pg_pz(3 * i + 0, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS - 6) = 1;
            pg_pz(3 * i + 0, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS + 0) = 1;
            pg_pz(3 * i + 1, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS - 5) = 1;
            pg_pz(3 * i + 1, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS + 1) = 1;
            pg_pz(3 * i + 2, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS - 4) = 1;
            pg_pz(3 * i + 2, lambdaOffset + i * (NUM_DEPENDENT_JOINTS + 6) + NUM_DEPENDENT_JOINTS + 2) = 1;
        }
    }
}

void CustomizedConstraints::compute_bounds() {
    g_lb = VecX::Zero(3 * N);
    g_ub = VecX::Zero(3 * N);
}

void CustomizedConstraints::print_violation_info() {
    for (int i = 0; i < N; i++) {
        const VecX& net_force = g.segment(3 * i, 3);
        if (net_force.norm() > 1e-3) {
            std::cout << "Violation at sample " << i << ": " << net_force.transpose() << std::endl;
        }
    }
}

}; // namespace DigitWholeBodySysID
}; // namespace RAPTOR
