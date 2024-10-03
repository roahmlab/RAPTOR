#include "BaseParametersIdentification.h"

namespace RAPTOR
{
namespace Kinova
{

    // // constructor
    // BaseParametersIdentification::BaseParametersIdentification()
    // {
    // }

    // // destructorS
    // BaseParametersIdentification::~BaseParametersIdentification()
    // {
    // }


    bool BaseParametersIdentification::set_parameters(
        MatX &Wh,
        VecX &Th,
        VecX &X,
        bool include_friction_offset,
        Model &model_input,
        std::shared_ptr<QRDecompositionSolver> regroupPtr,
        VecX &lb,
        VecX &ub,
        int b_full,
        int fm_dim,
        int Alg_case){ 

        Wh_=Wh;
        Th_=Th;
        X_=X;
        fm_dim_=fm_dim;
        include_friction_offset_=include_friction_offset;
        model_input_=model_input;
        lb_=lb;
        ub_=ub;
        b_full_ = b_full;
        Alg_case_ = Alg_case;
        nLinks_ = (Alg_case ==0)? model_input_.nv : 1;

        // regroup parameters 
        regroupPtr_ = regroupPtr;
        if(Alg_case ==0){
            Ginv_ = regroupPtr_->Ginv;
            pi_d_ = regroupPtr_->pi_d;
            b_dim_ = regroupPtr_->dim_id;
            d_dim_ = regroupPtr_->dim_d;
        }

        enable_hessian = false;
        return true;
    }

    bool BaseParametersIdentification::get_nlp_info(
        Index &n,
        Index &m,
        Index &nnz_jac_g,
        Index &nnz_h_lag,
        IndexStyleEnum &index_style)
    {
        // number of decision variables (independent and friction)
        n = b_full_;
        numVars= n;

        if (Alg_case_ !=2){
            m =nLinks_;
            numCons= m;
        }else{
            m =0;
            numCons= 0;
        }

        nnz_jac_g = n * m;
        nnz_h_lag = n * (n + 1) / 2;

        // use the C style indexing (0-based)
        index_style = TNLP::C_STYLE;

        return true;
    }

    bool BaseParametersIdentification::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                                   Index m, Number* g_l, Number* g_u)
    {
        // Set variable bounds
        for (Index i = 0; i < n; ++i)
        {
            x_l[i] = lb_[i];
            x_u[i] = ub_[i];
        }

        // Set constraint bounds
        for (Index i = 0; i < m; ++i)
        {
            g_l[i] = -1e19;
            g_u[i] = 0;
        }

        return true;
    }

    bool BaseParametersIdentification::get_starting_point(Index n, bool init_x, Number *x,
                                                                bool init_z, Number *z_L, Number *z_U,
                                                                Index m, bool init_lambda,
                                                                Number *lambda)
    {
        // Initialize x
        init_x = true;
        init_z = false;
        init_lambda = false;

        int interia_num = b_full_-fm_dim_*nLinks_;
        // Initialize intertia parameters after regroup
        for (Index i = 0; i < interia_num; ++i)
        {
            x[i] = X_[i];
        }

        // Initialize Ia
        Index idx = interia_num;
        for (Index i = 0; i< nLinks_; ++i){
            x[idx] = X_[idx];
        }
        
        // f_dim = fm_dim -Ia_dim
        for (Index j = 0; j < (fm_dim_-1) *nLinks_; ++j)
        {
            // x[idx] = X_[idx];
            x[idx] = 1;   /// change later 
            idx++;
        }

        return true;
    }


    bool BaseParametersIdentification::eval_f(
        Index n,
        const Number *x,
        bool new_x,
        Number &obj_value)
    {
        // Extract base parameters and dynamic parameters from x
        VecX z = Utils::initializeEigenVectorFromArray(x, n);

        VecX pi_s;
        if (Alg_case_ ==0){
            VecX pi_b = z.head(b_dim_);
            VecX pi_mf = z.tail(fm_dim_*nLinks_);
            pi_s.resize(b_dim_ + fm_dim_*nLinks_);
            pi_s << pi_b, pi_mf;
        }else{
            pi_s = z;
        }

        // Compute the objective function: ||Wh_ * [πb; πfm] - Th_||^2
        VecX diff = Wh_ * pi_s - Th_;
        obj_value = diff.squaredNorm();

        return true;
    }

    bool BaseParametersIdentification::eval_grad_f(
        Index n,
        const Number *x,
        bool new_x,
        Number *grad_f)
    {
        // Extract base parameters and dynamic parameters from x
        VecX z = Utils::initializeEigenVectorFromArray(x, n);
        VecX pi_s;
        if (Alg_case_ ==0){
            VecX pi_b = z.head(b_dim_);
            VecX pi_mf = z.tail(fm_dim_*nLinks_);
            pi_s.resize(b_dim_ + fm_dim_*nLinks_);
            pi_s << pi_b, pi_mf;
        }else{
            pi_s = z;
        }

        // // Grad: 2*W_h.T *(Wh_ * [πb; πfm] - Th_)
        VecX diff = Wh_ * pi_s - Th_;
        VecX grad = 2 * Wh_.transpose() * diff;

        // Assign gradient to grad_f
        for (Index i = 0; i < n; ++i)
        {
            grad_f[i] = grad[i];
        }

        return true;
    }

    bool BaseParametersIdentification::eval_g(Index n, const Number *x, bool new_x,
                                                Index m, Number *g)
    {
        if(Alg_case_ ==2){
            return true;
        }
        // Compute πs based on x
        VecX z = Utils::initializeEigenVectorFromArray(x, n);
        VecX pi_inertia;
        if (Alg_case_ ==0){
            VecX pi_b = z.head(b_dim_);
            VecX pi_mf = z.tail(fm_dim_*nLinks_);

            VecX pi_bd(b_dim_ + d_dim_);
            pi_bd << pi_b, pi_d_;
            pi_inertia = Ginv_ * pi_bd;
    
        }else{
            pi_inertia = z;
        }

        for (Index j = 0; j < nLinks_; ++j)
        {
            // Compute LMI matrix for link j
            MatX LMI(4,4);
            compute_LMI_matrix(pi_inertia, j, LMI);
            // Compute determinant
            double det = LMI.determinant();
            g[j] = -det;
        }

        return true;
    }

    bool BaseParametersIdentification::eval_jac_g(Index n, const Number *x, bool new_x,
                                                    Index m, Index nele_jac, Index *iRow, Index *jCol,
                                                    Number *values)
    {   
        if(Alg_case_ ==2){
            return true;
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
        }else {
            // Compute πs based on x
            VecX z = Utils::initializeEigenVectorFromArray(x, n);
            VecX pi_inertia;
            if (Alg_case_ ==0){
                VecX pi_b = z.head(b_dim_);
                VecX pi_mf = z.tail(fm_dim_*nLinks_);

                VecX pi_bd(b_dim_ + d_dim_);
                pi_bd << pi_b, pi_d_;
                pi_inertia = Ginv_ * pi_bd;
    
            }else{
                pi_inertia = z;
            }

            MatX DC;
            if (Alg_case_ ==0){
                DC.resize(nLinks_, 10*nLinks_);
            }else if (Alg_case_ ==1){
                DC.resize(nLinks_, (10+fm_dim_)*nLinks_);
            }
            DC.setZero();

            for (Index j = 0; j < nLinks_; ++j)
            {
                MatX LMI(4,4);
                compute_LMI_matrix(pi_inertia, j, LMI);
                double det = LMI.determinant();

                MatX dLMIdpi_full(4,4*10);// one joint
                compute_LMI_gradient(pi_inertia, j, dLMIdpi_full);

                 MatX LMI_inv;
                if (std::abs(det) < 1e-10) {
                    LMI_inv = LMI.completeOrthogonalDecomposition().pseudoInverse();
                } else {
                    LMI_inv = LMI.inverse();
                }

                for (int k = 0; k < 10; ++k) {
                    double grad_det = det * (LMI_inv * dLMIdpi_full.block<4,4>(0, k*4)).trace();
                    DC(j, j*10 + k) = -grad_det;  //  -det <= 0
                }

            }
            // calculate dphi_dX
            MatX dphi_dX;
            MatX s1(b_dim_ + d_dim_, b_dim_+ fm_dim_*nLinks_);
            // MatX L1;
            s1.setZero();
            s1.topLeftCorner(b_dim_, b_dim_).setIdentity();
            dphi_dX = Ginv_ * s1;
      
            MatX DC_regroup;
            if (Alg_case_ ==0){
                DC_regroup = DC * dphi_dX;
            }else{
                DC_regroup = DC;
            }
   
            // vector form, in row order
            // VecX grad;q
            Index idx = 0;
            for (Index row = 0; row < m; ++row) {
                for (Index col = 0; col < n; ++col) {
                    values[idx] = DC_regroup(row, col);
                    ++idx;
                }
            }
        }
        return true;
    }


    // bool BaseParametersIdentification::eval_h(
    //     Index n, const Number* x, bool new_x, Number obj_factor,
    //     Index m, const Number* lambda, bool new_lambda,
    //     Index nele_hess, Index* iRow, Index* jCol, Number* values)
    // {
    //     if (values == NULL) {
    //         return true;
    //     }
    //     for (Index idx = 0; idx < nele_hess; idx++) {
    //         values[idx] = 0.0;
    //     }

    //     return true;
    // }

    void BaseParametersIdentification::compute_LMI_matrix(const VecX &pi_inertia, Index j, MatX &LMI)
    {
        // Extract inertia tensor for link j
        Mat3 I_j;
        I_j << pi_inertia(10 * j), pi_inertia(10 * j + 1), pi_inertia(10 * j + 2),
               pi_inertia(10 * j + 1), pi_inertia(10 * j + 3), pi_inertia(10 * j + 4),
               pi_inertia(10 * j + 2), pi_inertia(10 * j + 4), pi_inertia(10 * j + 5);

        // center of mass
        Eigen::Vector3d l_j;
        l_j << pi_inertia(10 * j + 6), pi_inertia(10 * j + 7), pi_inertia(10 * j + 8);

        // mass
        Number M_j = pi_inertia(10 * j + 9);

        // Compute LMI matrix components
        Mat3 lmi11 = (I_j.trace() / 2.0) * Mat3::Identity() - I_j;
        Vec3 lmi12 = l_j;
        Number lmi22 = M_j;

        // Assemble LMI matrix
        LMI.setZero();
        LMI.block<3, 3>(0, 0) = lmi11;
        LMI.block<3, 1>(0, 3) = lmi12;
        LMI.block<1, 3>(3, 0) = lmi12.transpose();
        LMI(3, 3) = lmi22;
    }

    void BaseParametersIdentification::compute_LMI_gradient(const VecX &pi_full, Index j, MatX &dLMIdpi_full)
    {
        // Compute gradient of determinant with respect to pi_full
        MatX dLMIdpi_full1,
            dLMIdpi_full2,
            dLMIdpi_full3,
            dLMIdpi_full4,
            dLMIdpi_full5,
            dLMIdpi_full6,
            dLMIdpi_full7,
            dLMIdpi_full8,
            dLMIdpi_full9,
            dLMIdpi_full10;

        dLMIdpi_full1.resize(4, 4);
        dLMIdpi_full1.setZero();
        dLMIdpi_full1(0, 0) = -0.5;
        dLMIdpi_full1(1, 1) = 0.5;
        dLMIdpi_full1(2, 2) = 0.5;

        dLMIdpi_full2.resize(4, 4);
        dLMIdpi_full2.setZero();
        dLMIdpi_full2(0, 1) = -1;
        dLMIdpi_full2(1, 0) = -1;

        dLMIdpi_full3.resize(4, 4);
        dLMIdpi_full3.setZero();
        dLMIdpi_full3(0, 2) = -1;
        dLMIdpi_full3(2, 0) = -1;

        dLMIdpi_full4.resize(4, 4);
        dLMIdpi_full4.setZero();
        dLMIdpi_full4(0, 0) = 0.5;
        dLMIdpi_full4(1, 1) = -0.5;
        dLMIdpi_full4(2, 2) = 0.5;

        dLMIdpi_full5.resize(4, 4);
        dLMIdpi_full5.setZero();
        dLMIdpi_full5(1, 2) = -1;
        dLMIdpi_full5(2, 1) = -1;

        dLMIdpi_full6.resize(4, 4);
        dLMIdpi_full6.setZero();
        dLMIdpi_full6(0, 0) = 0.5;
        dLMIdpi_full6(1, 1) = 0.5;
        dLMIdpi_full6(2, 2) = -0.5;

        dLMIdpi_full7.resize(4, 4);
        dLMIdpi_full7.setZero();
        dLMIdpi_full7(0, 3) = 1;
        dLMIdpi_full7(3, 0) = 1;

        dLMIdpi_full8.resize(4, 4);
        dLMIdpi_full8.setZero();
        dLMIdpi_full8(1, 3) = 1;
        dLMIdpi_full8(3, 1) = 1;

        dLMIdpi_full9.resize(4, 4);
        dLMIdpi_full9.setZero();
        dLMIdpi_full9(2, 3) = 1;
        dLMIdpi_full9(3, 2) = 1;

        dLMIdpi_full10.resize(4, 4);
        dLMIdpi_full10.setZero();
        dLMIdpi_full10(3, 3) = 1;

        dLMIdpi_full << dLMIdpi_full1, dLMIdpi_full2, dLMIdpi_full3, dLMIdpi_full4, dLMIdpi_full5, dLMIdpi_full6, dLMIdpi_full7, dLMIdpi_full8, dLMIdpi_full9, dLMIdpi_full10;
    }


} // namespace Kinova
} // namespace RAPTOR
