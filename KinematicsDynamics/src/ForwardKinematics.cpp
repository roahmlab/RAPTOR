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

    T_identity = Transform();

    // allocate memory
    T = Transform();
    dTdq.resize(modelPtr_->nv);
    ddTddq.resize(modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        ddTddq[i].resize(modelPtr_->nv);
    }
    dddTdddq.resize(modelPtr_->nv);
    for (int i = 0; i < modelPtr_->nv; i++) {
        dddTdddq[i].resize(modelPtr_->nv);
        for (int j = 0; j < modelPtr_->nv; j++) {
            dddTdddq[i][j].resize(modelPtr_->nv);
        }
    }

    T_chains_collection.reserve(modelPtr_->nv * (modelPtr_->nv + 1) / 2);
    dTjdq_collections.reserve(modelPtr_->nv);
    ddTjddq_collections.reserve(modelPtr_->nv);
    dddTjdddq_collections.reserve(modelPtr_->nv);
}

const Transform& ForwardKinematicsSolver::getTransformDerivative(const VecX& q,
                                                                 const int id,
                                                                 const int order) {
    if (modelPtr_ == nullptr) {
        throw std::runtime_error("modelPtr_ is not initialized!");
    }

    if (jtype.size() != modelPtr_->nv) {
        throw std::runtime_error("jtype is not initialized!");
    }

    if (id < 0 || id >= modelPtr_->nv) {
        throw std::invalid_argument("id is out of range!");
    }

    if (order < 1 || order > 3) {
        throw std::invalid_argument("order has to be 1, 2, or 3!");
    }

    if (!Utils::ifTwoVectorEqual(current_q, q, 0)) {
        current_q = q;
        T_chains_collection.clear();
        dTjdq_collections.clear();
        ddTjddq_collections.clear();
        dddTjdddq_collections.clear();
    }

    const auto& jointPlacement = modelPtr_->jointPlacements[id + 1];

    if (order == 1) {
        if (dTjdq_collections.find(id) == dTjdq_collections.end()) {
            dTjdq_collections[id] = jointPlacement * Transform(jtype(id), q(id), 1);
        }

        return dTjdq_collections[id];
    }
    else if (order == 2) {
        if (ddTjddq_collections.find(id) == ddTjddq_collections.end()) {
            ddTjddq_collections[id] = jointPlacement * Transform(jtype(id), q(id), 2);
        }
        return ddTjddq_collections[id];
    }
    else if (order == 3) {
        if (dddTjdddq_collections.find(id) == dddTjdddq_collections.end()) {
            dddTjdddq_collections[id] = jointPlacement * Transform(jtype(id), q(id), 3);
        }
        return dddTjdddq_collections[id];
    }
    else {
        throw std::invalid_argument("order has to be less than or equal to 3!");
    }
}

const Transform& ForwardKinematicsSolver::getTransformChain(const VecX& q,
															const int chain_start,
															const int chain_end) {
    if (modelPtr_ == nullptr) {
        throw std::runtime_error("modelPtr_ is not initialized!");
    }

    if (jtype.size() != modelPtr_->nv) {
        throw std::runtime_error("jtype is not initialized!");
    }
    
    if (chain.empty()) {
        throw std::invalid_argument("chain is empty!");
    }

    if (chain_start > chain_end) {
        // throw std::invalid_argument("chain_start has to be less than chain_end!");
        return T_identity;
    }

    if (chain_start < 0 || chain_start >= chain.size() || 
        chain_end   < 0 || chain_end   >= chain.size()) {
        throw std::invalid_argument("chain_start or chain_end is out of range!");
    }

    const std::string chain_name = std::to_string(chain_start) + "-" + std::to_string(chain_end);

    if (!Utils::ifTwoVectorEqual(current_q, q, 0)) {
        current_q = q;
        T_chains_collection.clear();
        dTjdq_collections.clear();
        ddTjddq_collections.clear();
        dddTjdddq_collections.clear();
    }

    if (T_chains_collection.find(chain_name) != T_chains_collection.end()) {
        return T_chains_collection[chain_name];
    }

    if (chain_start == chain_end) {
        const int joint_id = chain[chain_start];
        const auto& jointPlacement = modelPtr_->jointPlacements[joint_id + 1];
        T_chains_collection[chain_name] = jointPlacement * 
                                          Transform(jtype(joint_id), q(joint_id));
        return T_chains_collection[chain_name];
    }

    // check the nearest chain
    const int joint_id = chain[chain_start];
    const auto& jointPlacement = modelPtr_->jointPlacements[joint_id + 1];
    T_chains_collection[chain_name] = jointPlacement * 
                                      Transform(jtype(joint_id), q(joint_id)) *
                                      getTransformChain(q, chain_start + 1, chain_end);
    
    return T_chains_collection[chain_name];
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

    if (startT != nullptr) {
        T_start = *startT;
    }
    else {
        T_start = Transform();
    }

    T = T_start * getTransformChain(q, 0, chain.size() - 1);

    if (order >= 1) {
        for (int i = 0; i < chain.size(); i++) {
            dTdq[chain[i]] = T_start *
                             getTransformChain(q, 0, i - 1) * 
                             getTransformDerivative(q, chain[i], 1) * 
                             getTransformChain(q, i + 1, chain.size() - 1);
        }
    }

    if (order >= 2) {
        for (int i = 0; i < chain.size(); i++) {
            for (int j = i; j < chain.size(); j++) {
                if (i == j) {
                    ddTddq[chain[i]][chain[j]] = T_start * 
                                                 getTransformChain(q, 0, i - 1) * 
                                                 getTransformDerivative(q, chain[i], 2) * 
                                                 getTransformChain(q, i + 1, chain.size() - 1);
                }
                else {
                    ddTddq[chain[i]][chain[j]] = T_start * 
                                                 getTransformChain(q, 0, i - 1) * 
                                                 getTransformDerivative(q, chain[i], 1) * 
                                                 getTransformChain(q, i + 1, j - 1) * 
                                                 getTransformDerivative(q, chain[j], 1) * 
                                                 getTransformChain(q, j + 1, chain.size() - 1);
                    ddTddq[chain[j]][chain[i]] = ddTddq[chain[i]][chain[j]];
                }
            }
        }
    }
    
    if (order >= 3) {
        for (int i = 0; i < chain.size(); i++) {
            for (int j = i; j < chain.size(); j++) {
                for (int k = j; k < chain.size(); k++) {
                    if (i == j && j == k) {
                        dddTdddq[chain[i]][chain[j]][chain[k]] = T_start * 
                                                                 getTransformChain(q, 0, i - 1) * 
                                                                 getTransformDerivative(q, chain[i], 3) * 
                                                                 getTransformChain(q, i + 1, chain.size() - 1);
                    }
                    else if (i == j) {
                        dddTdddq[chain[i]][chain[j]][chain[k]] = T_start * 
                                                                 getTransformChain(q, 0, i - 1) * 
                                                                 getTransformDerivative(q, chain[i], 2) *
                                                                 getTransformChain(q, i + 1, k - 1) * 
                                                                 getTransformDerivative(q, chain[k], 1) * 
                                                                 getTransformChain(q, k + 1, chain.size() - 1);
                        dddTdddq[chain[k]][chain[j]][chain[i]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[i]][chain[k]][chain[j]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                    }
                    else if (j == k) {
                        dddTdddq[chain[i]][chain[j]][chain[k]] = T_start * 
                                                                 getTransformChain(q, 0, i - 1) * 
                                                                 getTransformDerivative(q, chain[i], 1) * 
                                                                 getTransformChain(q, i + 1, j - 1) * 
                                                                 getTransformDerivative(q, chain[j], 2) * 
                                                                 getTransformChain(q, j + 1, chain.size() - 1);
                        dddTdddq[chain[j]][chain[i]][chain[k]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[j]][chain[k]][chain[i]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                    }
                    else {
                        dddTdddq[chain[i]][chain[j]][chain[k]] = T_start * 
                                                                 getTransformChain(q, 0, i - 1) * 
                                                                 getTransformDerivative(q, chain[i], 1) * 
                                                                 getTransformChain(q, i + 1, j - 1) * 
                                                                 getTransformDerivative(q, chain[j], 1) * 
                                                                 getTransformChain(q, j + 1, k - 1) * 
                                                                 getTransformDerivative(q, chain[k], 1) * 
                                                                 getTransformChain(q, k + 1, chain.size() - 1);
                        dddTdddq[chain[k]][chain[j]][chain[i]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[k]][chain[i]][chain[j]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[j]][chain[k]][chain[i]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[j]][chain[i]][chain[k]] = dddTdddq[chain[i]][chain[j]][chain[k]];
                        dddTdddq[chain[i]][chain[k]][chain[j]] = dddTdddq[chain[i]][chain[j]][chain[k]];
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