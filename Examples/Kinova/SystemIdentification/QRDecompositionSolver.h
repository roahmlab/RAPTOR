#ifndef QRDECOMPOSITIONSOLVER_H
#define QRDECOMPOSITIONSOLVER_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "RegressorInverseDynamics.h"
#include "Trajectories.h" 

namespace RAPTOR {
namespace Kinova {

class QRDecompositionSolver{
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    MatX Aid;
    MatX Ad;
    MatX Kd;
    MatX Ginv;
    VecX Beta;
    MatX RegroupMatrix;

    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
    // Eigen::VectorXi jtype = Eigen::VectorXi(0);


    int dim_id;
    int dim_d;
    double  eps = 1e-8;

    QRDecompositionSolver();
    QRDecompositionSolver(const MatX& ObservationMatrix, const VecX& InParam);
    ~QRDecompositionSolver();

private:
    void getdata();
    void compute(const MatX& W, const Eigen::VectorXd& InParam);
};

}; // namespace Kinova
}; // namespace RAPTOR



#endif