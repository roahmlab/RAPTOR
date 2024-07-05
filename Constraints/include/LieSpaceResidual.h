#ifndef LIE_SPACE_RESIDUAL_H    
#define LIE_SPACE_RESIDUAL_H    

#include <unsupported/Eigen/MatrixFunctions>

#include "Utils.h"
#include "ForwardKinematics.h"

namespace IDTO {

// This namespace contains functions to compute the residual of the Lie space
// including the translation and rotation part in SE(3) space, with their derivatives
// Assume that fkPtr_ is a valid pointer and has computed all corresponding kinematics
namespace LieSpaceResidual {

double safedacosdx(const double x,
                   const bool throw_exception = false);

double safeddacosddx(const double x,
                     const bool throw_exception = false);

// compute x / sin(x) safely and accurately around x = 0
double safexSinx(const double x,
                 const double nearZeroThreshold = 1e-6);

// compute d(x / sin(x)) / dx (the derivative of the previous function x / sin(x)) 
// safely and accurately around x = 0
double safedxSinxdx(const double x,
                    const double nearZeroThreshold = 1e-6);

// compute dd(x / sin(x)) / ddx (the derivative of the previous function x / sin(x)) 
// safely and accurately around x = 0
double safeddxSinxddx(const double x,
                      const double nearZeroThreshold = 1e-6);

// You have to make sure that fkPtr_ has already computed the corresponding forward kinematics
// and its gradient or hessian before using this function!!!
Eigen::Vector3d translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3d& desiredPosition,
                                    Eigen::MatrixXd* gradientPtr_ = nullptr,
                                    Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_ = nullptr);

// You have to make sure that fkPtr_ has already computed the corresponding forward kinematics
// and its gradient or hessian before using this function!!!
Eigen::Vector3d rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3d& desiredRotation,
                                 Eigen::MatrixXd* gradientPtr_ = nullptr,
                                 Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_ = nullptr);

}; // namespace LieSpaceResidual
}; // namespace IDTO

#endif // LIE_SPACE_RESIDUAL_H