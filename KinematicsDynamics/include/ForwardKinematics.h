#ifndef FORWARD_KINEMATICS_HPP
#define FORWARD_KINEMATICS_HPP

#include "Utils.h"
#include "Transform.h"
#include "HigherOrderDerivatives.h"

namespace RAPTOR {

class ForwardKinematicsSolver {
public:
    using Model = pinocchio::Model;
	using Vec3 = Eigen::Vector3d;
	using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
	using MatX = Eigen::MatrixXd;

    // Constructor
	ForwardKinematicsSolver() = default;

	ForwardKinematicsSolver(const Model* model_input,
							Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
	~ForwardKinematicsSolver() = default;

    // class methods:
		// compute forward kinematics and its high-order derivatives
	void compute(const int start,
				 const int end,
				 const VecX& q, 
				 const Transform* startT = nullptr,
				 const Transform* endT = nullptr,
				 const int order = 0);

	// get the forward kinematics result in different formats
	Transform getTransform() const;

	const Transform& getTransformChain(const VecX& q,
									   const int start,
									   const int end);

	const Transform& getTransformDerivative(const VecX& q,
											const int id,
											const int order);

	Vec3 getTranslation() const;

	Mat3 getRotation() const;

	Vec3 getRPY() const; // roll, pitch, yaw

	MatX getTranslationJacobian() const;

	void getRotationJacobian(Eigen::Array<Mat3, Eigen::Dynamic, 1>& result) const;

	MatX getRPYJacobian() const;

	void getTranslationHessian(Eigen::Array<MatX, 3, 1>& result) const;

	void getRotationHessian(Eigen::Array<MatX, 3, 3>& result) const;
	void getRotationHessian(Eigen::Array<Mat3, Eigen::Dynamic, Eigen::Dynamic>& result) const;

	void getRPYHessian(Eigen::Array<MatX, 3, 1>& result) const;

	// return third-order-tensor * x, which is a 2D matrix
	void getTranslationThirdOrderTensor(const VecX& x, Eigen::Array<MatX, 3, 1>& result) const;

	// return third-order-tensor * x, which is a 2D matrix
	void getRPYThirdOrderTensor(const VecX& x, Eigen::Array<MatX, 3, 1>& result) const;
	
    // class members:
	const Model* modelPtr_ = nullptr;
	Eigen::VectorXi jtype;

		// a index vector that stores the kinematics chain
	std::vector<int> chain;  

		// results
	Transform T;
	std::vector<Transform> dTdq;
	std::vector<std::vector<Transform>> ddTddq;
	std::vector<std::vector<std::vector<Transform>>> dddTdddq;

	Transform T_identity;
	Transform T_start;

	VecX current_q;
	std::unordered_map<std::string, Transform> T_chains_collection;
	std::unordered_map<int, Transform> dTjdq_collections;
	std::unordered_map<int, Transform> ddTjddq_collections;
	std::unordered_map<int, Transform> dddTjdddq_collections;
};

}  // namespace RAPTOR

#endif