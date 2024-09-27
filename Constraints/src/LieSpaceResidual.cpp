#include "LieSpaceResidual.h"

namespace RAPTOR {
namespace LieSpaceResidual {

Eigen::Vector3d translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3d& desiredPosition,
                                    Eigen::MatrixXd* gradientPtr_,
                                    Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_) {
    if (gradientPtr_ != nullptr) {
        *gradientPtr_ = fkPtr_->getTranslationJacobian();
    }

    if (hessianPtr_ != nullptr) {
        fkPtr_->getTranslationHessian(*hessianPtr_);
    }
    
    return fkPtr_->getTranslation() - desiredPosition;
}

Eigen::Vector3d rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3d& desiredRotation,
                                 Eigen::MatrixXd* gradientPtr_,
                                 Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_) {
    bool compute_derivatives = (gradientPtr_ != nullptr) ||
                               (hessianPtr_ != nullptr);
    bool compute_hessian = (hessianPtr_ != nullptr);

    // kinematics chain (derivative is only related to these joints)
    const auto& chain = fkPtr_->chain;

    const Eigen::Matrix3d currentRotation = fkPtr_->getRotation();
    Eigen::Matrix3d residualMatrix = desiredRotation.transpose() * currentRotation;

    Eigen::Array<Eigen::Matrix3d, Eigen::Dynamic, 1> dRdq;
    Eigen::Array<Eigen::Matrix3d, Eigen::Dynamic, Eigen::Dynamic> ddRddq;
    Eigen::Tensor<Eigen::Matrix3d, 3> dddRdddq;
    if (compute_derivatives) {
        fkPtr_->getRotationJacobian(dRdq);

        for (auto i : chain) {
            dRdq(i) = desiredRotation.transpose() * dRdq(i);
        }

        if (compute_hessian) {
            fkPtr_->getRotationHessian(ddRddq);

            for (auto i : chain) {
                for (auto j : chain) {
                    ddRddq(i, j) = desiredRotation.transpose() * ddRddq(i, j);
                }
            }
        }
    }

    double traceR = residualMatrix.trace();

    Eigen::VectorXd dtraceRdq;
    Eigen::MatrixXd ddtraceRddq;
    Eigen::Tensor<double, 3> dddtraceRdddq;
    if (compute_derivatives) {
        dtraceRdq = Eigen::VectorXd::Zero(dRdq.size());
        for (auto i : chain) {
            dtraceRdq(i) = dRdq(i).trace();
        }

        if (compute_hessian) {
            ddtraceRddq = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddtraceRddq(i, j) = ddRddq(i, j).trace();
                }
            }
        }
    }

    double theta = HigherOrderDerivatives::safeacos((traceR - 1) / 2);

    Eigen::VectorXd dthetadq;
    Eigen::MatrixXd ddthetaddq;
    Eigen::Tensor<double, 3> dddthetadddq;
    if (compute_derivatives) {
        const double dacosdxTraceR = HigherOrderDerivatives::safedacosdx((traceR - 1) / 2);
        
        dthetadq = Eigen::VectorXd::Zero(dRdq.size());
        for (auto i : chain) {
            dthetadq(i) = 0.5 * dacosdxTraceR * dtraceRdq(i);
        }

        if (compute_hessian) {
            const double ddacosddxTraceR = HigherOrderDerivatives::safeddacosddx((traceR - 1) / 2);

            ddthetaddq = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddthetaddq(i, j) = 0.5 * 
                        (0.5 * ddacosddxTraceR * dtraceRdq(i) * dtraceRdq(j) +
                         dacosdxTraceR * ddtraceRddq(i, j));
                }
            }
        }
    }

    const double st = sin(theta);
    const double ct = cos(theta);
    const double xSinxTheta = HigherOrderDerivatives::safexSinx(theta);
    const Eigen::Matrix3d RRT = residualMatrix - residualMatrix.transpose();

    Eigen::Matrix3d logR = 0.5 * xSinxTheta * RRT;

    if (compute_derivatives) {
        const double dxSinxTheta = HigherOrderDerivatives::safedxSinxdx(theta);

        if (gradientPtr_ != nullptr) {
            Eigen::MatrixXd& gradient = *gradientPtr_;
            gradient = Eigen::MatrixXd::Zero(3, dRdq.size());

            for (auto i : chain) {
                Eigen::Matrix3d temp1 = 
                    dxSinxTheta * dthetadq(i) * RRT;
                Eigen::Matrix3d temp2 = 
                    xSinxTheta * (dRdq(i) - dRdq(i).transpose());
                gradient.col(i) = Utils::unskew(0.5 * (temp1 + temp2));
            }
        }

        if (compute_hessian) {
            const double ddxSinxddx = HigherOrderDerivatives::safeddxSinxddx(theta);

            Eigen::Array<Eigen::MatrixXd, 3, 1>& hessian = *hessianPtr_;
            hessian(0) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(1) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(2) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());

            for (auto i : chain) {
                for (auto j : chain) {
                    Eigen::Matrix3d temp1_1 = 
                        ddxSinxddx * dthetadq(i) * dthetadq(j) * RRT;
                    Eigen::Matrix3d temp1_2 = 
                        dxSinxTheta * ddthetaddq(i, j) * RRT;
                    Eigen::Matrix3d temp2 = 
                        dxSinxTheta * dthetadq(i) * (dRdq(j) - dRdq(j).transpose());
                    Eigen::Matrix3d temp3 = 
                        dxSinxTheta * dthetadq(j) * (dRdq(i) - dRdq(i).transpose());
                    Eigen::Matrix3d temp4 = 
                        xSinxTheta * (ddRddq(i, j) - ddRddq(i, j).transpose());

                    Eigen::Vector3d h = 
                        Utils::unskew(
                            0.5 * ((temp1_1 + temp1_2) + temp2 + temp3 + temp4));

                    hessian(0)(i, j) = h(0);
                    hessian(1)(i, j) = h(1);
                    hessian(2)(i, j) = h(2);
                }
            }
        }
    }
    
    return Utils::unskew(logR);
}

}; // namespace LieSpaceResidual
}; // namespace RAPTOR