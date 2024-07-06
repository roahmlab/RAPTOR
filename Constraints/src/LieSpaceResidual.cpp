#include "LieSpaceResidual.h"

namespace IDTO {
namespace LieSpaceResidual {

Eigen::Vector3d translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3d& desiredPosition,
                                    Eigen::MatrixXd* gradientPtr_,
                                    Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_,
                                    Eigen::Array<Eigen::Tensor<double, 3>, 3, 1>* thridOrderTensorPtr_) {
    if (gradientPtr_ != nullptr) {
        *gradientPtr_ = fkPtr_->getTranslationJacobian();
    }

    if (hessianPtr_ != nullptr) {
        fkPtr_->getTranslationHessian(*hessianPtr_);
    }

    if (thridOrderTensorPtr_ != nullptr) {
        fkPtr_->getTranslationThirdOrderTensor(*thridOrderTensorPtr_);
    }
    
    return fkPtr_->getTranslation() - desiredPosition;
}

Eigen::Vector3d rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3d& desiredRotation,
                                 Eigen::MatrixXd* gradientPtr_,
                                 Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_,
                                 Eigen::Array<Eigen::Tensor<double, 3>, 3, 1>* thridOrderTensorPtr_) {
    bool compute_derivatives = (gradientPtr_ != nullptr) ||
                               (hessianPtr_ != nullptr) ||
                               (thridOrderTensorPtr_ != nullptr);
    bool compute_hessian = (hessianPtr_ != nullptr) ||
                           (thridOrderTensorPtr_ != nullptr);
    bool compute_thrid_order_tensor = (thridOrderTensorPtr_ != nullptr);

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

            if (compute_thrid_order_tensor) {
                fkPtr_->getRotationThirdOrderTensor(dddRdddq);

                for (auto i : chain) {
                    for (auto j : chain) {
                        for (auto k : chain) {
                            dddRdddq(i, j, k) = desiredRotation.transpose() * dddRdddq(i, j, k);
                        }
                    }
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

            if (compute_thrid_order_tensor) {
                dddtraceRdddq.resize(dRdq.size(), dRdq.size(), dRdq.size());
                dddtraceRdddq.setZero();
                for (auto i : chain) {
                    for (auto j : chain) {
                        for (auto k : chain) {
                            dddtraceRdddq(i, j, k) = dddRdddq(i, j, k).trace();
                        }
                    }
                }
            }
        }
    }

    double theta = Utils::safeacos((traceR - 1) / 2);

    Eigen::VectorXd dthetadq;
    Eigen::MatrixXd ddthetaddq;
    Eigen::Tensor<double, 3> dddthetadddq;
    if (compute_derivatives) {
        const double dacosdxTraceR = safedacosdx((traceR - 1) / 2);
        
        dthetadq = Eigen::VectorXd::Zero(dRdq.size());
        for (auto i : chain) {
            dthetadq(i) = 0.5 * dacosdxTraceR * dtraceRdq(i);
        }

        if (compute_hessian) {
            const double ddacosddxTraceR = safeddacosddx((traceR - 1) / 2);

            ddthetaddq = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddthetaddq(i, j) = 0.5 * 
                        (0.5 * ddacosddxTraceR * dtraceRdq(i) * dtraceRdq(j) +
                         dacosdxTraceR * ddtraceRddq(i, j));
                }
            }

            if (compute_thrid_order_tensor) {
                const double dddacosdddxTraceR = safedddacosdddx((traceR - 1) / 2);

                dddthetadddq.resize(dRdq.size(), dRdq.size(), dRdq.size());
                dddthetadddq.setZero();
                for (auto i : chain) {
                    for (auto j : chain) {
                        for (auto k : chain) {
                            dddthetadddq(i, j, k) = 0.5 * 
                                (0.25 * dddacosdddxTraceR * dtraceRdq(i) * dtraceRdq(j) * dtraceRdq(k) +
                                 0.5 * ddacosddxTraceR * ddtraceRddq(i, k) * dtraceRdq(j) +
                                 0.5 * ddacosddxTraceR * dtraceRdq(i) * ddtraceRddq(j, k) +
                                 0.5 * ddacosddxTraceR * dtraceRdq(k) * ddtraceRddq(i, j) +
                                 dacosdxTraceR * dddtraceRdddq(i, j, k));
                        }
                    }
                }
            }
        }
    }

    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const double xSinxTheta = safexSinx(theta);
    const Eigen::Matrix3d RRT = residualMatrix - residualMatrix.transpose();

    Eigen::Matrix3d logR = 0.5 * xSinxTheta * RRT;

    if (compute_derivatives) {
        const double dxSinxTheta = safedxSinxdx(theta);

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
            const double ddxSinxddx = safeddxSinxddx(theta);

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

double safedacosdx(const double x,
                   const bool throw_exception) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return -1e10; // a very large negative number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    return -1.0 / std::sqrt(1.0 - x * x);
}

double safeddacosddx(const double x,
                     const bool throw_exception) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return -1e10; // a very large negative number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    return -x / std::pow(1.0 - x * x, 1.5);
}

double safedddacosdddx(const double x,
                       const bool throw_exception) {
    if (x >= 1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is greater than 1.0");
        }
        return -1e10; // a very large negative number
    } 
    else if (x <= -1.0) {
        if (throw_exception) {
            throw std::runtime_error("Input value is less than -1.0");
        }
        return -1e10; // a very large negative number
    } 
    const double xSquare = x * x;
    return -(1 + 2 * xSquare) / std::pow(1.0 - xSquare, 2.5);
}

double safexSinx(const double x,
                 const double nearZeroThreshold) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xFourth = xSquare * xSquare;
        return 1.0 + xSquare / 6.0 + 7.0 * xFourth / 360.0; // + O(x^6)
    }
    return x / std::sin(x);
}

double safedxSinxdx(const double x,
                    double nearZeroThreshold) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xThird = x * xSquare;
        double XFifth = xThird * xSquare;
        return x / 3.0 + 7.0 * xThird / 90.0 + 31.0 * XFifth / 2520.0; // + O(x^7)
    }
    const double sinx = std::sin(x);
    return (sinx - x * std::cos(x)) / (sinx * sinx);
}

double safeddxSinxddx(const double x,
                      const double nearZeroThreshold) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xFourth = xSquare * xSquare;
        double xSixth = xFourth * xSquare;
        return 1.0 / 3.0 + 7.0 * xSquare / 30.0 + 31.0 * xFourth / 504.0; // + O(x^6)
    }
    const double sinx = std::sin(x);
    const double cosx = std::cos(x);
    const double sinxSquare = sinx * sinx;
    return (x * sinxSquare - 2 * cosx * sinx + 2 * x * cosx * cosx) / (sinxSquare * sinx);
}

double safedddxSinxdddx(const double x,
                        const double nearZeroThreshold) {
    if (fabs(x) < nearZeroThreshold) { // use Taylor expansion to approximate
        double xSquare = x * x;
        double xThird = x * xSquare;
        double XFifth = xThird * xSquare;
        return 7.0 * x / 15.0 + 31.0 * xThird / 126.0 + 127.0 * XFifth / 1800.0; // + O(x^7)
    }
    const double sinx = std::sin(x);
    const double cosx = std::cos(x);
    const double sinxSquare = sinx * sinx;
    return x * cosx / sinxSquare +
           6 * cosx * cosx / (sinx * sinxSquare) -
           6 * x * cosx / (sinxSquare * sinxSquare) + 
           3 / sinx;
}

}; // namespace LieSpaceResidual
}; // namespace IDTO