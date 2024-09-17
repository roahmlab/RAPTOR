#ifndef LIE_SPACE_RESIDUAL_H    
#define LIE_SPACE_RESIDUAL_H    

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Utils.h"
#include "HigherOrderDerivatives.h"
#include "ForwardKinematics.h"

namespace RAPTOR {

// This namespace contains functions to compute the residual of the Lie space
// including the translation and rotation part in SE(3) space, with their derivatives
// Assume that fkPtr_ is a valid pointer and has computed all corresponding kinematics
namespace LieSpaceResidual {

// You have to make sure that fkPtr_ has already computed the corresponding forward kinematics
// and its gradient or hessian before using this function!!!
Eigen::Vector3f translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3f& desiredPosition,
                                    Eigen::MatrixXf* gradientPtr_ = nullptr,
                                    Eigen::Array<Eigen::MatrixXf, 3, 1>* hessianPtr_ = nullptr);

// You have to make sure that fkPtr_ has already computed the corresponding forward kinematics
// and its gradient or hessian before using this function!!!
Eigen::Vector3f rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3f& desiredRotation,
                                 Eigen::MatrixXf* gradientPtr_ = nullptr,
                                 Eigen::Array<Eigen::MatrixXf, 3, 1>* hessianPtr_ = nullptr);

}; // namespace LieSpaceResidual
}; // namespace RAPTOR

#endif // LIE_SPACE_RESIDUAL_H