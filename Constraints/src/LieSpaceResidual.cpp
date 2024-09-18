#include "LieSpaceResidual.h"

namespace RAPTOR {
namespace LieSpaceResidual {

Eigen::Vector3f translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3f& desiredPosition,
                                    Eigen::MatrixXf* gradientPtr_,
                                    Eigen::Array<Eigen::MatrixXf, 3, 1>* hessianPtr_) {
    if (gradientPtr_ != nullptr) {
        *gradientPtr_ = fkPtr_->getTranslationJacobian();
    }

    if (hessianPtr_ != nullptr) {
        fkPtr_->getTranslationHessian(*hessianPtr_);
    }
    
    return fkPtr_->getTranslation() - desiredPosition;
}

Eigen::Vector3f rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3f& desiredRotation,
                                 Eigen::MatrixXf* gradientPtr_,
                                 Eigen::Array<Eigen::MatrixXf, 3, 1>* hessianPtr_) {
    bool compute_derivatives = (gradientPtr_ != nullptr) ||
                               (hessianPtr_ != nullptr);
    bool compute_hessian = (hessianPtr_ != nullptr);

    // kinematics chain (derivative is only related to these joints)
    const auto& chain = fkPtr_->chain;

    const Eigen::Matrix3f currentRotation = fkPtr_->getRotation();
    Eigen::Matrix3f residualMatrix = desiredRotation.transpose() * currentRotation;

    Eigen::Array<Eigen::Matrix3f, Eigen::Dynamic, 1> dRdq;
    Eigen::Array<Eigen::Matrix3f, Eigen::Dynamic, Eigen::Dynamic> ddRddq;
    Eigen::Tensor<Eigen::Matrix3f, 3> dddRdddq;
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

    float traceR = residualMatrix.trace();

    Eigen::VectorXf dtraceRdq;
    Eigen::MatrixXf ddtraceRddq;
    Eigen::Tensor<float, 3> dddtraceRdddq;
    if (compute_derivatives) {
        dtraceRdq = Eigen::VectorXf::Zero(dRdq.size());
        for (auto i : chain) {
            dtraceRdq(i) = dRdq(i).trace();
        }

        if (compute_hessian) {
            ddtraceRddq = Eigen::MatrixXf::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddtraceRddq(i, j) = ddRddq(i, j).trace();
                }
            }
        }
    }

    float theta = HigherOrderDerivatives::safeacos((traceR - 1) / 2);

    Eigen::VectorXf dthetadq;
    Eigen::MatrixXf ddthetaddq;
    Eigen::Tensor<float, 3> dddthetadddq;
    if (compute_derivatives) {
        const float dacosdxTraceR = HigherOrderDerivatives::safedacosdx((traceR - 1) / 2);
        
        dthetadq = Eigen::VectorXf::Zero(dRdq.size());
        for (auto i : chain) {
            dthetadq(i) = 0.5 * dacosdxTraceR * dtraceRdq(i);
        }

        if (compute_hessian) {
            const float ddacosddxTraceR = HigherOrderDerivatives::safeddacosddx((traceR - 1) / 2);

            ddthetaddq = Eigen::MatrixXf::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddthetaddq(i, j) = 0.5 * 
                        (0.5 * ddacosddxTraceR * dtraceRdq(i) * dtraceRdq(j) +
                         dacosdxTraceR * ddtraceRddq(i, j));
                }
            }
        }
    }

    const float st = sinf(theta);
    const float ct = cosf(theta);
    const float xSinxTheta = HigherOrderDerivatives::safexSinx(theta);
    const Eigen::Matrix3f RRT = residualMatrix - residualMatrix.transpose();

    Eigen::Matrix3f logR = 0.5 * xSinxTheta * RRT;

    if (compute_derivatives) {
        const float dxSinxTheta = HigherOrderDerivatives::safedxSinxdx(theta);

        if (gradientPtr_ != nullptr) {
            Eigen::MatrixXf& gradient = *gradientPtr_;
            gradient = Eigen::MatrixXf::Zero(3, dRdq.size());

            for (auto i : chain) {
                Eigen::Matrix3f temp1 = 
                    dxSinxTheta * dthetadq(i) * RRT;
                Eigen::Matrix3f temp2 = 
                    xSinxTheta * (dRdq(i) - dRdq(i).transpose());
                gradient.col(i) = Utils::unskew(0.5 * (temp1 + temp2));
            }
        }

        if (compute_hessian) {
            const float ddxSinxddx = HigherOrderDerivatives::safeddxSinxddx(theta);

            Eigen::Array<Eigen::MatrixXf, 3, 1>& hessian = *hessianPtr_;
            hessian(0) = Eigen::MatrixXf::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(1) = Eigen::MatrixXf::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(2) = Eigen::MatrixXf::Zero(ddRddq.rows(), ddRddq.cols());

            for (auto i : chain) {
                for (auto j : chain) {
                    Eigen::Matrix3f temp1_1 = 
                        ddxSinxddx * dthetadq(i) * dthetadq(j) * RRT;
                    Eigen::Matrix3f temp1_2 = 
                        dxSinxTheta * ddthetaddq(i, j) * RRT;
                    Eigen::Matrix3f temp2 = 
                        dxSinxTheta * dthetadq(i) * (dRdq(j) - dRdq(j).transpose());
                    Eigen::Matrix3f temp3 = 
                        dxSinxTheta * dthetadq(j) * (dRdq(i) - dRdq(i).transpose());
                    Eigen::Matrix3f temp4 = 
                        xSinxTheta * (ddRddq(i, j) - ddRddq(i, j).transpose());

                    Eigen::Vector3f h = 
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