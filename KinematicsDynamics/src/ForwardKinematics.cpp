#include "ForwardKinematics.h"

namespace IDTO {

// ForwardKinematicsSolver::ForwardKinematicsSolver() {
//     q_copy = VecX::Zero(0);
// }

ForwardKinematicsSolver::ForwardKinematicsSolver(const Model* model_input,
										         const Eigen::VectorXi& jtype_input) : 
    modelPtr_(model_input),
    jtype(jtype_input) {
    if (modelPtr_->nv != jtype.size()) {
        std::cerr << "modelPtr_->nv = " << modelPtr_->nv << std::endl;
        std::cerr << "jtype.size() = " << jtype.size() << std::endl;
        throw std::invalid_argument("modelPtr_->nv != jtype.size()");
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

        if (search_id < 0 || search_id > modelPtr_->nv) {
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

}  // namespace IDTO