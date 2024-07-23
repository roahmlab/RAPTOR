#include "Optimizer.h"
#include "ForwardKinematics.h"

using namespace RAPTOR;
using namespace Ipopt;

class FKGradientChecker : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    FKGradientChecker() = default;

    /** Default destructor */
    ~FKGradientChecker() = default;

    bool set_parameters(
        const VecX& x0_input,
        const Model& model_input, 
        const Eigen::VectorXi& jtype_input
    ) {
        enable_hessian = true;
        x0 = x0_input;
        model = model_input;
        jtype = jtype_input;

        fkPtr_ = std::make_unique<ForwardKinematicsSolver>(&model, jtype);

        return true;
    }

    bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    ) final override
    {
        // number of decision variables
        numVars = model.nq;
        n = numVars;

        m = 0;

        nnz_jac_g = n * m;
        nnz_h_lag = n * (n + 1) / 2;

        // use the C style indexing (0-based)
        index_style = TNLP::C_STYLE;

        return true;
    }

    bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    ) final override 
    {
        if(n != numVars){
            THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_f!");
        }

        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        // fkPtr_->compute(0, model.nq, z);
        // Vec3 rpy = fkPtr_->getRPY();

        // obj_value = roll_weight * rpy(0) + 
        //             pitch_weight * rpy(1) + 
        //             yaw_weight * rpy(2);

        fkPtr_->compute(0, model.nq, z, nullptr, nullptr, 1);
        MatX Jrpy = fkPtr_->getRPYJacobian();
        Vec3 rpy = Jrpy * VecX::Ones(n);

        obj_value = roll_weight * rpy(0) + 
                    pitch_weight * rpy(1) + 
                    yaw_weight * rpy(2);

        return true;
    }

    bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    ) final override
    {
        if(n != numVars){
            THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_grad_f!");
        }

        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        // fkPtr_->compute(0, model.nq, z, nullptr, nullptr, 1);
        // MatX Jrpy = fkPtr_->getRPYJacobian();

        // for ( Index i = 0; i < n; i++ ) {
        //     grad_f[i] = roll_weight * Jrpy(0, i) + 
        //                 pitch_weight * Jrpy(1, i) + 
        //                 yaw_weight * Jrpy(2, i);
        // }

        fkPtr_->compute(0, model.nq, z, nullptr, nullptr, 2);
        Eigen::Array<MatX, 3, 1> Hrpy;
        fkPtr_->getRPYHessian(Hrpy);

        VecX Jr = Hrpy(0) * VecX::Ones(n);
        VecX Jp = Hrpy(1) * VecX::Ones(n);
        VecX Jy = Hrpy(2) * VecX::Ones(n);

        for ( Index i = 0; i < n; i++ ) {
            grad_f[i] = roll_weight * Jr(i) + 
                        pitch_weight * Jp(i) + 
                        yaw_weight * Jy(i);
        }

        return true;
    }

    bool eval_hess_f(
        Index         n,
        const Number* x,
        bool          new_x,
        MatX&         hess_f
    ) final override
    {
        if(n != numVars){
            THROW_EXCEPTION(IpoptException, "*** Error wrong value of n in eval_hess_f!");
        }

        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        // fkPtr_->compute(0, model.nq, z, nullptr, nullptr, 2);
        // Eigen::Array<MatX, 3, 1> Hrpy;
        // fkPtr_->getRPYHessian(Hrpy);

        // hess_f = roll_weight * Hrpy(0) + 
        //          pitch_weight * Hrpy(1) + 
        //          yaw_weight * Hrpy(2);

        fkPtr_->compute(0, model.nq, z, nullptr, nullptr, 3);
        Eigen::Array<MatX, 3, 1> Trpy;
        fkPtr_->getRPYThirdOrderTensor(VecX::Ones(n), Trpy);

        hess_f = roll_weight * Trpy(0) + 
                 pitch_weight * Trpy(1) + 
                 yaw_weight * Trpy(2);

        return true;
    }
    
    FKGradientChecker(
       const FKGradientChecker&
    );

    FKGradientChecker& operator=(
       const FKGradientChecker&
    );

    Model model;
    Eigen::VectorXi jtype;

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_;

    const double roll_weight = 1.0;
    const double pitch_weight = 1.0;
    const double yaw_weight = 1.0;
};

int main() {
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    Eigen::VectorXd z0 = Eigen::VectorXd::Random(model.nq);

    // Initialize gradient checker
    SmartPtr<FKGradientChecker> mynlp = new FKGradientChecker();
    try {
	    mynlp->set_parameters(z0,
                              model,
                              jtype);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	app->Options()->SetNumericValue("max_wall_time", 1e-5);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("hessian_approximation", "exact");

    // For gradient checking
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("derivative_test", "second-order");
    app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    try {
        status = app->OptimizeTNLP(mynlp);
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}