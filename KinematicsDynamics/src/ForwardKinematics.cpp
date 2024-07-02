#include "ForwardKinematics.h"

namespace IDTO {

// ForwardKinematicsSolver::ForwardKinematicsSolver() {
//     q_copy = VecX::Zero(0);
// }

ForwardKinematicsSolver::ForwardKinematicsSolver(const Model& model_input,
										         const Eigen::VectorXi& jtype_input) : 
    model(model_input),
    jtype(jtype_input) {
    if (model.nv != jtype.size()) {
        std::cerr << "model.nv = " << model.nv << std::endl;
        std::cerr << "jtype.size() = " << jtype.size() << std::endl;
        throw std::invalid_argument("model.nv != jtype.size()");
    }
}

void ForwardKinematicsSolver::compute(const int start,
                                      const int end,
                                      const VecX& q, 
                                      const Transform* startT,
                                      const Transform* endT,
                                      const int order) {
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

        if (search_id < 0 || search_id > model.nv) {
            throw std::runtime_error("Can not find the end joint in the model!");
        }

        search_id = model.parents[search_id];
    }
    std::reverse(chain.begin(), chain.end());

    // allocate memory
    T = Transform();
    if (order >= 1) {
        dTdq.resize(model.nv);

        if (order >= 2) {
            ddTddq.resize(model.nv);
            for (int i = 0; i < model.nv; i++) {
                ddTddq[i].resize(model.nv);
            }

            if (order >= 3) {
                dddTdddq.resize(model.nv);
                for (int i = 0; i < model.nv; i++) {
                    dddTdddq[i].resize(model.nv);
                    for (int j = 0; j < model.nv; j++) {
                        dddTdddq[i][j].resize(model.nv);
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

    auto start_clock = std::chrono::high_resolution_clock::now();
    // iterative process to compute the forward kinematics
    for (auto i : chain) {
        // pinocchio joint index starts from 1
        const auto& jointPlacement = model.jointPlacements[i + 1];
        
        Transform Tj(jtype(i), q(i));
        // Transform dTjdq(jtype(i), q(i), 1);
        // Transform ddTjddq(jtype(i), q(i), 2);
        // Transform dddTjdddq(jtype(i), q(i), 3);

        T *= (jointPlacement * Tj);
        
        // if (order >= 1) {
        //     for (auto j : chain) {
        //         dTdq[j] *= jointPlacement;
        //         if (j == i) {
        //             dTdq[j] *= dTjdq;
        //         }
        //         else {
        //             dTdq[j] *= Tj;
        //         }

        //         if (order >= 2) {
        //             for (auto k : chain) {
        //                 if (k >= j) {
        //                     ddTddq[j][k] *= jointPlacement;
        //                     if (j == i && k == i) {
        //                         ddTddq[j][k] *= ddTjddq;
        //                     } 
        //                     else if (j == i || k == i) {
        //                         ddTddq[j][k] *= dTjdq;
        //                     } 
        //                     else {
        //                         ddTddq[j][k] *= Tj;
        //                     }
        //                 } 
        //                 else {
        //                     ddTddq[j][k] = ddTddq[k][j];
        //                 }

        //                 if (order >= 3) {
        //                     for (auto h : chain) {
        //                         if (h >= k && k >= j) {
        //                             dddTdddq[j][k][h] *= jointPlacement;
        //                             if (j == i && k == i && h == i) {
        //                                 dddTdddq[j][k][h] *= dddTjdddq;
        //                             } 
        //                             else if ((j == i && k == i) || (j == i && h == i) || (k == i && h == i)) {
        //                                 dddTdddq[j][k][h] *= ddTjddq;
        //                             } 
        //                             else if (j == i || k == i || h == i) {
        //                                 dddTdddq[j][k][h] *= dTjdq;
        //                             }
        //                             else {
        //                                 dddTdddq[j][k][h] *= Tj;
        //                             }
        //                         }
        //                         else {
        //                             if (h < k) {
        //                                 dddTdddq[j][k][h] = dddTdddq[j][h][k];
        //                             }
        //                             else if (k < j) {
        //                                 dddTdddq[j][k][h] = dddTdddq[k][j][h];
        //                             }
        //                             else {
        //                                 dddTdddq[j][k][h] = dddTdddq[h][k][j];
        //                             }
        //                         }
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }
    }
    auto stop_clock = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_clock - start_clock);
    std::cout << "Internal FK: " << duration.count() << " nanoseconds" << std::endl;

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

Eigen::Vector3d ForwardKinematicsSolver::getTranslation() const {
    return Eigen::Vector3d(T.data[9], 
                           T.data[10], 
                           T.data[11]);
}

Eigen::Matrix3d ForwardKinematicsSolver::getRotation() const {
    Eigen::Matrix3d R;
    R << T.data[0], T.data[1], T.data[2],
         T.data[3], T.data[4], T.data[5],
         T.data[6], T.data[7], T.data[8];
    return R;
}

Eigen::MatrixXd ForwardKinematicsSolver::getTranslationJacobian() const {
    if (dTdq.size() != model.nv) {
        throw std::runtime_error("dTdq is not computed yet!");
    }

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, model.nv);
    for (auto i : chain) {
        J.col(i) = Eigen::Vector3d(dTdq[i].data[9], 
                                   dTdq[i].data[10], 
                                   dTdq[i].data[11]);
    }
    return J;
}

void ForwardKinematicsSolver::getTranslationHessian(Eigen::Array<MatX, 3, 1>& result) const {
    if (ddTddq.size() != model.nv) {
        throw std::runtime_error("ddTddq is not computed yet!");
    }

    for (int i = 0; i < 3; i++) {
        result(i) = Eigen::MatrixXd::Zero(model.nv, model.nv);
    }

    for (auto i : chain) {
        for (auto j : chain) {
            result(0)(i, j) = ddTddq[i][j].data[9];
            result(1)(i, j) = ddTddq[i][j].data[10];
            result(2)(i, j) = ddTddq[i][j].data[11];
        }
    }
}

}  // namespace IDTO