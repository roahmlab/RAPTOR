#ifndef QRDECOMPOSITIONSOLVER_H
#define QRDECOMPOSITIONSOLVER_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "RegressorInverseDynamics.h"
#include "Trajectories.h"
#include "Polynomials.h"

namespace RAPTOR {

class QRDecompositionSolver {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Output matrices and vectors
    MatX Aid;
    MatX Ad;
    MatX A;
    MatX Kd;
    MatX Ginv;
    VecX Beta;
    MatX RegroupMatrix;
    VecX pi_id;
    VecX pi_d;
    VecX InputParams;

    // Smart pointers to Trajectories and RegressorInverseDynamics
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
    Eigen::VectorXi jtype;

    // Dimensions of independent and dependent parameters
    int dim_id;
    int dim_d;

    // Threshold for numerical stability
    double eps;

    // Constructors and Destructor
    QRDecompositionSolver();
    QRDecompositionSolver(const MatX& ObservationMatrix, const VecX& InputParams);
    ~QRDecompositionSolver();

    // Member functions
    void getData();
    void compute(const MatX& W, const VecX& InputParams);
};

}; // namespace RAPTOR

#endif // QRDECOMPOSITIONSOLVER_H
