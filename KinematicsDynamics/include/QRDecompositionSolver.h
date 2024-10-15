#ifndef QR_DECOMPOSITION_SOLVER_H
#define QR_DECOMPOSITION_SOLVER_H

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "pinocchio/algorithm/regressor.hpp"

namespace RAPTOR {

class QRDecompositionSolver {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;

    // Constructors 
    QRDecompositionSolver() = default;

    QRDecompositionSolver(const Model& model_input);

    // Destructor
    ~QRDecompositionSolver() = default;

    // class methods:
        // generate random observation on a certain number of instances
    void generateRandomObservation(const int numInstances = 1000);

        // compute QR decomposition-based regrouping
    void computeRegroupMatrix(const double eps = 1e-8);

    // class members:
    std::shared_ptr<Model> modelPtr_ = nullptr;
    std::shared_ptr<Data> dataPtr_ = nullptr;

    VecX phi; // original inertial parameters (ungrouped)
    MatX ObservationMatrix; // observation matrix

        // intermediate matrices and vectors
    MatX Aid;
    MatX Ad;
    MatX A;
    MatX Kd;
    MatX Ginv;
    VecX phi_id; // independent parameters
    VecX phi_d; // dependent parameters

        // results are stored in these matrices and vectors
    VecX beta; // regrouped inertial parameters
    MatX RegroupMatrix;

    // Dimensions of independent and dependent parameters
    int dim_id = 0;
    int dim_d = 0;
};

}; // namespace RAPTOR

#endif // QR_DECOMPOSITION_SOLVER_H
