#include "Optimizer.h"

namespace IDTO {

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool Optimizer::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
    if (n != numVars) {
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
    }

    // fill in a Eigen Vector instance of decision variables
    VecX z(n);
    for ( Index i = 0; i < n; i++ ) {
        z(i) = x[i];
    }

    Index iter = 0;
    for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
        // compute constraints
        try {
            constraintsPtrVec_[c]->compute(z, false);
        }
        catch (std::exception& e) {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("*** Error in eval_g!");
        }

        // fill in constraints
        for ( Index i = 0; i < constraintsPtrVec_[c]->m; i++ ) {
            g[iter] = constraintsPtrVec_[c]->g(i);
            iter++;
        }
    }

    return true;
}
// [TNLP_eval_g]


// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool Optimizer::eval_jac_g(
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
        throw std::runtime_error("*** Error wrong value of n in eval_g!");
    }
    if (m != numCons) {
        throw std::runtime_error("*** Error wrong value of m in eval_g!");
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
        // fill in a Eigen Vector instance of decision variables
        VecX z(n);
        for ( Index i = 0; i < n; i++ ) {
            z(i) = x[i];
        }

        Index iter = 0;
        for (Index c = 0; c < constraintsPtrVec_.size(); c++) {
            // compute constraints
            try {
                constraintsPtrVec_[c]->compute(z, true);
            }
            catch (std::exception& e) {
                std::cout << e.what() << std::endl;
                throw std::runtime_error("*** Error in eval_jac_g!");
            }

            // fill in constraints
            for ( Index i = 0; i < constraintsPtrVec_[c]->pg_pz.rows(); i++ ) {
                for ( Index j = 0; j < constraintsPtrVec_[c]->pg_pz.cols(); j++ ) {
                    values[iter] = constraintsPtrVec_[c]->pg_pz(i, j);
                    iter++;
                }
            }
        }
    }

    return true;
}
// [TNLP_eval_jac_g]


// [TNLP_eval_h]
//return the structure or values of the Hessian
bool Optimizer::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
    return false;
}
// [TNLP_eval_h]

}; // namespace IDTO