#ifndef FORWARD_KINEMATICS_HPP
#define FORWARD_KINEMATICS_HPP

#include <unsupported/Eigen/CXX11/Tensor>
#include "Transform.h"
#include "HigherOrderDerivatives.h"

namespace IDTO {

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
							const Eigen::VectorXi& jtype_input);

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

	void getTranslationThirdOrderTensor(Eigen::Array<Eigen::Tensor<double, 3>, 3, 1>& result) const;

	void getRotationThirdOrderTensor(Eigen::Tensor<Mat3, 3>& result) const;
	
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

	// 	// internal copies
	// int current_order = -1;
	// int end_copy = -1;
	// int start_copy = -1;
	// VecX q_copy;
	// Transform T_copy;
	// std::vector<Transform> dTdq_copy;
	// std::vector<std::vector<Transform>> ddTddq_copy;
	// std::vector<std::vector<std::vector<Transform>>> dddTdddq_copy;
};

}  // namespace IDTO

#endif