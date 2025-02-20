#include "ForwardKinematics.h"

namespace RAPTOR {

ForwardKinematicsSolver::ForwardKinematicsSolver(const Model* model_input,
										         Eigen::VectorXi jtype_input) : 
    modelPtr_(model_input),
    jtype(jtype_input) {
    if (jtype.size() > 0) {
        if (modelPtr_->nv != jtype.size()) {
            std::cerr << "modelPtr_->nv = " << modelPtr_->nv << std::endl;
            std::cerr << "jtype.size() = " << jtype.size() << std::endl;
            throw std::invalid_argument("modelPtr_->nv != jtype.size()");
        }
    }
    else {
        jtype = convertPinocchioJointType(*modelPtr_);
    }
}

void ForwardKinematicsSolver::compute(const int start,
                                      const int end,
                                      const VecX& q, 
                                      const Transform* startT,
                                      const Transform* endT,
                                      const int order) {
    if (modelPtr_ == nullptr) {
        throw std::runtime_error("modelPtr_ is not initialized!");
    }

    if (jtype.size() != modelPtr_->nv) {
        throw std::runtime_error("jtype is not initialized!");
    }

    if (order < 0) {
        throw std::invalid_argument("order has to be non-negative!");
    }
    if (order > 3) {
        throw std::invalid_argument("order has to be less than or equal to 3!");
    }

    // find the kinematics chain
    chain.clear();
    int search_id = end;
    while (search_id != start) {
        // pinocchio joint index starts from 1
        chain.push_back(search_id - 1);

        if (search_id < 0 || search_id > modelPtr_->nv || chain.size() > modelPtr_->nv) {
            throw std::runtime_error("Can not find the end joint in the modelPtr_!");
        }

        search_id = modelPtr_->parents[search_id];
    }
    std::reverse(chain.begin(), chain.end());

    // allocate memory
    T = Transform();
    if (order >= 1) {
        dTdq.resize(modelPtr_->nv);

        if (order >= 2) {
            ddTddq.resize(modelPtr_->nv);
            for (int i = 0; i < modelPtr_->nv; i++) {
                ddTddq[i].resize(modelPtr_->nv);
            }

            if (order >= 3) {
                dddTdddq.resize(modelPtr_->nv);
                for (int i = 0; i < modelPtr_->nv; i++) {
                    dddTdddq[i].resize(modelPtr_->nv);
                    for (int j = 0; j < modelPtr_->nv; j++) {
                        dddTdddq[i][j].resize(modelPtr_->nv);
                    }
                }
            }
        }
    }

    // initialize memory with start transformation matrix
    if (startT != nullptr) {
        T = (*startT);

        if (order >= 1) {
            for (auto i : chain) {
                dTdq[i] = (*startT);

                if (order >= 2) {
                    for (auto j : chain) {
                        ddTddq[i][j] = (*startT);

                        if (order >= 3) {
                            for (auto k : chain) {
                                dddTdddq[i][j][k] = (*startT);
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        T = Transform();

        if (order >= 1) {
            for (auto i : chain) {
                dTdq[i] = Transform();

                if (order >= 2) {
                    for (auto j : chain) {
                        ddTddq[i][j] = Transform();

                        if (order >= 3) {
                            for (auto k : chain) {
                                dddTdddq[i][j][k] = Transform();
                            }
                        }
                    }
                }
            }
        }
    }

    // iterative process to compute the forward kinematics
    for (auto i : chain) {
        // pinocchio joint index starts from 1
        const auto& jointPlacement = modelPtr_->jointPlacements[i + 1];
        
        Transform Tj(jtype(i), q(i));
        Transform dTjdq(jtype(i), q(i), 1);
        Transform ddTjddq(jtype(i), q(i), 2);
        Transform dddTjdddq(jtype(i), q(i), 3);

        T *= (jointPlacement * Tj);
        
        if (order >= 1) {
            for (auto j : chain) {
                dTdq[j] *= jointPlacement;
                if (j == i) {
                    dTdq[j] *= dTjdq;
                }
                else {
                    dTdq[j] *= Tj;
                }

                if (order >= 2) {
                    for (auto k : chain) {
                        if (k >= j) {
                            ddTddq[j][k] *= jointPlacement;
                            if (j == i && k == i) {
                                ddTddq[j][k] *= ddTjddq;
                            } 
                            else if (j == i || k == i) {
                                ddTddq[j][k] *= dTjdq;
                            } 
                            else {
                                ddTddq[j][k] *= Tj;
                            }
                        } 
                        else {
                            ddTddq[j][k] = ddTddq[k][j];
                        }

                        if (order >= 3) {
                            for (auto h : chain) {
                                if (h >= k && k >= j) {
                                    dddTdddq[j][k][h] *= jointPlacement;
                                    if (j == i && k == i && h == i) {
                                        dddTdddq[j][k][h] *= dddTjdddq;
                                    } 
                                    else if ((j == i && k == i) || (j == i && h == i) || (k == i && h == i)) {
                                        dddTdddq[j][k][h] *= ddTjddq;
                                    } 
                                    else if (j == i || k == i || h == i) {
                                        dddTdddq[j][k][h] *= dTjdq;
                                    }
                                    else {
                                        dddTdddq[j][k][h] *= Tj;
                                    }
                                }
                                else {
                                    if (h < k) {
                                        dddTdddq[j][k][h] = dddTdddq[j][h][k];
                                    }
                                    else if (k < j) {
                                        dddTdddq[j][k][h] = dddTdddq[k][j][h];
                                    }
                                    else {
                                        dddTdddq[j][k][h] = dddTdddq[h][k][j];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // multiply the end transformation matrix
    if (endT != nullptr) {
        T *= (*endT);

        if (order >= 1) {
            for (auto i : chain) {
                dTdq[i] *= (*endT);

                if (order >= 2) {
                    for (auto j : chain) {
                        ddTddq[i][j] *= (*endT);

                        if (order >= 3) {
                            for (auto k : chain) {
                                dddTdddq[i][j][k] *= (*endT);
                            }
                        }
                    }
                }
            }
        
        }
    }
}

Transform ForwardKinematicsSolver::getTransform() const {
    return T;
}

Eigen::Vector3d ForwardKinematicsSolver::getTranslation() const {
    return T.p;
}

Eigen::Matrix3d ForwardKinematicsSolver::getRotation() const {
    return T.R;
}

Eigen::Vector3d ForwardKinematicsSolver::getRPY() const {
    return T.getRPY();
}

Eigen::MatrixXd ForwardKinematicsSolver::getTranslationJacobian() const {
    if (dTdq.size() != modelPtr_->nv) {
        throw std::runtime_error("dTdq is not computed yet!");
    }

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, modelPtr_->nv);
    for (auto i : chain) {
        J.col(i) = dTdq[i].p;
    }
    
    return J;
}

void ForwardKinematicsSolver::getRotationJacobian(Eigen::Array<Mat3, Eigen::Dynamic, 1>& result) const {
    if (dTdq.size() != modelPtr_->nv) {
        throw std::runtime_error("dTdq is not computed yet!");
    }
    
    result.resize(modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        result(i).setZero();
    }

    for (auto i : chain) {
        result(i) = dTdq[i].R;
    }
}

Eigen::MatrixXd ForwardKinematicsSolver::getRPYJacobian() const {
    if (dTdq.size() != modelPtr_->nv) {
        throw std::runtime_error("dTdq is not computed yet!");
    }
    
    const Mat3& R = T.R; // so that the code is cleaner

    const double rollDenom = R(1,2) * R(1,2) + R(2,2) * R(2,2);
    const double yawDenom = R(0,0) * R(0,0) + R(0,1) * R(0,1);

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, modelPtr_->nv);
    for (auto i : chain) {
        const Mat3& dRdq = dTdq[i].R; // so that the code is cleaner

        J(0, i) = (R(1,2) * dRdq(2,2) - 
                   R(2,2) * dRdq(1,2)) 
                        / rollDenom;

        J(1, i) = HigherOrderDerivatives::safedasindx(R(0,2)) * dRdq(0,2);

        J(2, i) = (R(0,1) * dRdq(0,0) - 
                   R(0,0) * dRdq(0,1)) 
                        / yawDenom;
    }
    
    return J;
}

void ForwardKinematicsSolver::getTranslationHessian(Eigen::Array<MatX, 3, 1>& result) const {
    if (ddTddq.size() != modelPtr_->nv) {
        throw std::runtime_error("ddTddq is not computed yet!");
    }

    for (int i = 0; i < 3; i++) {
        result(i) = Eigen::MatrixXd::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    for (auto i : chain) {
        for (auto j : chain) {
            result(0)(i, j) = ddTddq[i][j].p(0);
            result(1)(i, j) = ddTddq[i][j].p(1);
            result(2)(i, j) = ddTddq[i][j].p(2);
        }
    }
}

void ForwardKinematicsSolver::getRotationHessian(Eigen::Array<MatX, 3, 3>& result) const {
    if (ddTddq.size() != modelPtr_->nv) {
        throw std::runtime_error("ddTddq is not computed yet!");
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result(i, j) = Eigen::MatrixXd::Zero(modelPtr_->nv, modelPtr_->nv);
        }
    }

    for (auto i : chain) {
        for (auto j : chain) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    result(k, l)(i, j) = ddTddq[i][j].R(k, l);
                }
            }
        }
    }
}

void ForwardKinematicsSolver::getRotationHessian(Eigen::Array<Mat3, Eigen::Dynamic, Eigen::Dynamic>& result) const {
    if (ddTddq.size() != modelPtr_->nv) {
        throw std::runtime_error("ddTddq is not computed yet!");
    }
    
    result.resize(modelPtr_->nv, modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        for (int j = 0; j < modelPtr_->nv; j++) {
            result(i, j).setZero();
        }
    }

    for (auto i : chain) {
        for (auto j : chain) {
            result(i, j) = ddTddq[i][j].R;
        }
    }
}

void ForwardKinematicsSolver::getRPYHessian(Eigen::Array<MatX, 3, 1>& result) const {
    if (ddTddq.size() != modelPtr_->nv) {
        throw std::runtime_error("ddTddq is not computed yet!");
    }

    for (int i = 0; i < 3; i++) {
        result(i) = Eigen::MatrixXd::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    const Mat3& R = T.R; // so that the code is cleaner

    const double R12Square = R(1,2) * R(1,2);
    const double R22Square = R(2,2) * R(2,2);
    const double rollDenom = R12Square + R22Square;
    const double rollDenomSquare = rollDenom * rollDenom;

    const double R00Square = R(0,0) * R(0,0);
    const double R01Square = R(0,1) * R(0,1);
    const double yawDenom =  R00Square + R01Square;
    const double yawDenomSquare = yawDenom * yawDenom;

    const double dasindx = HigherOrderDerivatives::safedasindx(R(0,2));
    const double ddasinddx = HigherOrderDerivatives::safeddasinddx(R(0,2));

    for (auto i : chain) {
        const Mat3& dRdqi = dTdq[i].R; // so that the code is cleaner
        for (auto j : chain) {
            const Mat3& dRdqj = dTdq[j].R; // so that the code is cleaner
            const Mat3& ddRdqdq = ddTddq[i][j].R; // so that the code is cleaner

            result(0)(i, j) = 
                dRdqj(1,2) * (dRdqi(2,2) / rollDenom + 
                              (R(1,2) * R(2,2) * dRdqi(1,2) - dRdqi(2,2) * R12Square)
                                 * (2.0 / rollDenomSquare)) -
                dRdqj(2,2) * (dRdqi(1,2) / rollDenom + 
                              (R(1,2) * R(2,2) * dRdqi(2,2) - dRdqi(1,2) * R22Square)
                                 * (2.0 / rollDenomSquare)) +
                (R(1,2) * ddRdqdq(2,2) - R(2,2) * ddRdqdq(1,2)) / rollDenom;

            result(1)(i, j) = 
                ddasinddx * dRdqi(0,2) * dRdqj(0,2) + 
                dasindx   * ddRdqdq(0,2);

            result(2)(i, j) = 
                dRdqj(0,1) * (dRdqi(0,0) / yawDenom + 
                              (R(0,1) * R(0,0) * dRdqi(0,1) - dRdqi(0,0) * R01Square)
                                 * (2.0 / yawDenomSquare)) -
                dRdqj(0,0) * (dRdqi(0,1) / yawDenom + 
                              (R(0,1) * R(0,0) * dRdqi(0,0) - dRdqi(0,1) * R00Square)
                                 * (2.0 / yawDenomSquare)) +
                (R(0,1) * ddRdqdq(0,0) - R(0,0) * ddRdqdq(0,1)) / yawDenom;
        }
    }
}

void ForwardKinematicsSolver::getTranslationThirdOrderTensor(const VecX& x, Eigen::Array<MatX, 3, 1>& result) const {
    if (dddTdddq.size() != modelPtr_->nv) {
        throw std::runtime_error("dddTdddq is not computed yet!");
    }

    if (x.size() != modelPtr_->nv) {
        throw std::invalid_argument("x has to be of size modelPtr_->nv");
    }

    for (int i = 0; i < 3; i++) {
        result(i) = Eigen::MatrixXd::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    for (auto i : chain) {
        for (auto j : chain) {
            for (auto k : chain) {
                result(0)(i, j) += dddTdddq[i][j][k].p(0) * x(k);
                result(1)(i, j) += dddTdddq[i][j][k].p(1) * x(k);
                result(2)(i, j) += dddTdddq[i][j][k].p(2) * x(k);
            }
        }
    }
}

void ForwardKinematicsSolver::getRPYThirdOrderTensor(const VecX& x, Eigen::Array<MatX, 3, 1>& result) const {
    if (dddTdddq.size() != modelPtr_->nv) {
        throw std::runtime_error("dddTdddq is not computed yet!");
    }

    if (x.size() != modelPtr_->nv) {
        throw std::invalid_argument("x has to be of size modelPtr_->nv");
    }

    for (int i = 0; i < 3; i++) {
        result(i) = Eigen::MatrixXd::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    const Mat3& R = T.R; // so that the code is cleaner

    const double R12Square = R(1,2) * R(1,2);
    const double R22Square = R(2,2) * R(2,2);
    const double rollDenom = R12Square + R22Square;
    const double rollDenomSquare = rollDenom * rollDenom;
    const double rollDenomThird = rollDenomSquare * rollDenom;

    const double R00Square = R(0,0) * R(0,0);
    const double R01Square = R(0,1) * R(0,1);
    const double yawDenom =  R00Square + R01Square;
    const double yawDenomSquare = yawDenom * yawDenom;
    const double yawDenomThird = yawDenomSquare * yawDenom;

    const double dasindx = HigherOrderDerivatives::safedasindx(R(0,2));
    const double ddasinddx = HigherOrderDerivatives::safeddasinddx(R(0,2));
    const double dddasindddx = HigherOrderDerivatives::safedddasindddx(R(0,2));

    for (auto i : chain) {
        const Mat3& dRdqi = dTdq[i].R; // so that the code is cleaner

        const double temp1 = 2.0 * R(1,2) * dRdqi(1,2) / rollDenomSquare;
        const double temp2 = 2.0 * R(2,2) * dRdqi(2,2) / rollDenomSquare;
        const double temp3 = 8.0 * R(1,2) * dRdqi(1,2) * R22Square / rollDenomThird;
        const double temp4 = 8.0 * R(2,2) * dRdqi(2,2) * R12Square / rollDenomThird;
        const double temp5 = temp1 - temp2 + temp4 - temp3;

        const double temp6 = 2.0 * R(0,1) * dRdqi(0,1) / yawDenomSquare;
        const double temp7 = 2.0 * R(0,0) * dRdqi(0,0) / yawDenomSquare;
        const double temp8 = 8.0 * R(0,1) * dRdqi(0,1) * R00Square / yawDenomThird;
        const double temp9 = 8.0 * R(0,0) * dRdqi(0,0) * R01Square / yawDenomThird;
        const double temp10 = temp6 - temp7 + temp9 - temp8;

        for (auto j : chain) {
            const Mat3& dRdqj = dTdq[j].R; // so that the code is cleaner
            const Mat3& ddRdqidqj = ddTddq[i][j].R; // so that the code is cleaner
            for (auto k : chain) {
                const Mat3& dRdqk = dTdq[k].R; // so that the code is cleaner
                const Mat3& ddRdqidqk = ddTddq[i][k].R; // so that the code is cleaner
                const Mat3& ddRdqjdqk = ddTddq[j][k].R; // so that the code is cleaner
                const Mat3& dddRdqdqdq = dddTdddq[i][j][k].R; // so that the code is cleaner

                const double T0_i_j_k = 
                    -ddRdqjdqk(2,2) * (R(1,2) * temp2 + dRdqi(1,2) / rollDenom - dRdqi(1,2) * R22Square / rollDenomSquare * 2.0) 
                    + ddRdqjdqk(1,2) * (R(2,2) * temp1 + dRdqi(2,2) / rollDenom - dRdqi(2,2) * R12Square / rollDenomSquare * 2.0)
                    + dRdqk(1,2) * (
                        -dRdqj(1,2) * (
                            R(1,2) * dRdqi(2,2) / rollDenomSquare * 6.0 
                            - R(2,2) * dRdqi(1,2) / rollDenomSquare * 2.0 
                            - (R(1,2) * R(1,2) * R(1,2)) * dRdqi(2,2) / rollDenomThird * 8.0 
                            + R(2,2) * dRdqi(1,2) * R12Square / rollDenomThird * 8.0
                        ) 
                        + dRdqj(2,2) * temp5 
                        + ddRdqidqj(2,2) / rollDenom 
                        - ddRdqidqj(2,2) * R12Square / rollDenomSquare * 2.0 
                        + R(1,2) * R(2,2) * ddRdqidqj(1,2) / rollDenomSquare * 2.0
                    ) 
                    - dRdqk(2,2) * (
                        dRdqj(2,2) * (
                            R(1,2) * dRdqi(2,2) / rollDenomSquare * 2.0 
                            - R(2,2) * dRdqi(1,2) / rollDenomSquare * 6.0 
                            + (R(2,2) * R(2,2) * R(2,2)) * dRdqi(1,2) / rollDenomThird * 8.0 
                            - R(1,2) * dRdqi(2,2) * R22Square / rollDenomThird * 8.0
                        ) 
                        - dRdqj(1,2) * temp5 
                        + ddRdqidqj(1,2) / rollDenom 
                        - ddRdqidqj(1,2) * R22Square / rollDenomSquare * 2.0 
                        + R(1,2) * R(2,2) * ddRdqidqj(2,2) / rollDenomSquare * 2.0
                    ) 
                    - ddRdqidqk(1,2) * (
                        dRdqj(2,2) * (1.0 / rollDenom - R22Square / rollDenomSquare * 2.0) 
                        - R(1,2) * R(2,2) * dRdqj(1,2) / rollDenomSquare * 2.0
                    ) 
                    + ddRdqidqk(2,2) * (
                        dRdqj(1,2) * (1.0 / rollDenom - R12Square / rollDenomSquare * 2.0) 
                        - R(1,2) * R(2,2) * dRdqj(2,2) / rollDenomSquare * 2.0
                    ) 
                    + R(1,2) * dddRdqdqdq(2,2) / rollDenom 
                    - R(2,2) * dddRdqdqdq(1,2) / rollDenom;


                const double T1_i_j_k = 
                    dddasindddx * dRdqi(0,2) * dRdqj(0,2) * dRdqk(0,2) +
                    ddasinddx * ddRdqidqk(0,2) * dRdqj(0,2) +
                    ddasinddx * ddRdqjdqk(0,2) * dRdqi(0,2) +
                    ddasinddx * ddRdqidqj(0,2) * dRdqk(0,2) +
                    dasindx * dddRdqdqdq(0,2);

                const double T2_i_j_k = 
                    -ddRdqjdqk(0,0) * (R(0,1) * temp7 + dRdqi(0,1) / yawDenom - dRdqi(0,1) * R00Square / yawDenomSquare * 2.0) 
                    + ddRdqjdqk(0,1) * (R(0,0) * temp6 + dRdqi(0,0) / yawDenom - dRdqi(0,0) * R01Square / yawDenomSquare * 2.0)
                    + dRdqk(0,1) * (
                        -dRdqj(0,1) * (
                            R(0,1) * dRdqi(0,0) / yawDenomSquare * 6.0 
                            - R(0,0) * dRdqi(0,1) / yawDenomSquare * 2.0 
                            - (R(0,1) * R(0,1) * R(0,1)) * dRdqi(0,0) / yawDenomThird * 8.0 
                            + R(0,0) * dRdqi(0,1) * R01Square / yawDenomThird * 8.0
                        ) 
                        + dRdqj(0,0) * temp10 
                        + ddRdqidqj(0,0) / yawDenom 
                        - ddRdqidqj(0,0) * R01Square / yawDenomSquare * 2.0 
                        + R(0,1) * R(0,0) * ddRdqidqj(0,1) / yawDenomSquare * 2.0
                    ) 
                    - dRdqk(0,0) * (
                        dRdqj(0,0) * (
                            R(0,1) * dRdqi(0,0) / yawDenomSquare * 2.0 
                            - R(0,0) * dRdqi(0,1) / yawDenomSquare * 6.0 
                            + (R(0,0) * R(0,0) * R(0,0)) * dRdqi(0,1) / yawDenomThird * 8.0 
                            - R(0,1) * dRdqi(0,0) * R00Square / yawDenomThird * 8.0
                        ) 
                        - dRdqj(0,1) * temp10 
                        + ddRdqidqj(0,1) / yawDenom 
                        - ddRdqidqj(0,1) * R00Square / yawDenomSquare * 2.0 
                        + R(0,1) * R(0,0) * ddRdqidqj(0,0) / yawDenomSquare * 2.0
                    ) 
                    - ddRdqidqk(0,1) * (
                        dRdqj(0,0) * (1.0 / yawDenom - R00Square / yawDenomSquare * 2.0) 
                        - R(0,1) * R(0,0) * dRdqj(0,1) / yawDenomSquare * 2.0
                    ) 
                    + ddRdqidqk(0,0) * (
                        dRdqj(0,1) * (1.0 / yawDenom - R01Square / yawDenomSquare * 2.0) 
                        - R(0,1) * R(0,0) * dRdqj(0,0) / yawDenomSquare * 2.0
                    ) 
                    + R(0,1) * dddRdqdqdq(0,0) / yawDenom 
                    - R(0,0) * dddRdqdqdq(0,1) / yawDenom;

                result(0)(i, j) += T0_i_j_k * x(k);
                result(1)(i, j) += T1_i_j_k * x(k);
                result(2)(i, j) += T2_i_j_k * x(k);
            }
        }
    }
}

}  // namespace RAPTOR