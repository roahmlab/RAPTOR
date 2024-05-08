#ifndef FORWARD_KINEMATICS_HPP
#define FORWARD_KINEMATICS_HPP

#include "Transform.h"

namespace IDTO {

class ForwardKinematicsHighOrderDerivative {
public:
    using Model = pinocchio::Model;
	using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;

    // Constructor
	ForwardKinematicsHighOrderDerivative();

    // Destructor
	~ForwardKinematicsHighOrderDerivative() = default;

    // class methods:
	void fk(Transform& T, 
			const Model& model, 
			const Eigen::VectorXi& jtype, 
			const int end, 
			const int start, 
			const VecX& q, 
			const Transform& endT, 
			const Transform& startT);

	void fk_jacobian(std::vector<Transform>& dTdq, 
					 const Model& model,
					 const Eigen::VectorXi& jtype,
					 const int end, 
					 const int start, 
					 const VecX& q, 
					 const Transform& endT, 
					 const Transform& startT);

	void fk_hessian(std::vector<std::vector<Transform>>& ddTddq, 
					const Model& model, 
					const Eigen::VectorXi& jtype,
					const int end,
					const int start, 
					const VecX& q, 
					const Transform& endT, 
					const Transform& startT);

	void fk_thirdorder(std::vector<std::vector<std::vector<Transform>>>& dddTdddq, 
					   const Model& model, 
					   const Eigen::VectorXi& jtype,
					   const int end,
					   const int start, 
					   const VecX& q, 
					   const Transform& endT, 
					   const Transform& startT);

	Vec3 Transform2xyz(const Transform& T);

	VecX Transform2xyzrpy(const Transform& T);

	void Transform2xyzJacobian(Eigen::MatrixXd& J, 
							   const Transform& T, 
							   const std::vector<Transform>& dTdq);

	void Transform2xyzrpyJacobian(Eigen::MatrixXd& J, 
								  const Transform& T, 
								  const std::vector<Transform>& dTdq);

    void Transform2xyzHessian(Eigen::Array<Eigen::MatrixXd, 3, 1>& H,
							  const Transform& T, 
							  const std::vector<Transform>& dTdq,
							  const std::vector<std::vector<Transform>>& ddTddq);

	void Transform2xyzrpyHessian(Eigen::Array<Eigen::MatrixXd, 6, 1>& H,
								 const Transform& T, 
								 const std::vector<Transform>& dTdq,
								 const std::vector<std::vector<Transform>>& ddTddq);

	void Transform2xyzrpyThirdOrder(Eigen::Array<Eigen::MatrixXd, 6, 1>& TOx,
									const VecX& x,
									const Transform& T, 
									const std::vector<Transform>& dTdq,
									const std::vector<std::vector<Transform>>& ddTddq,
									const std::vector<std::vector<std::vector<Transform>>>& dddTdddq);

    // class members:
    // a index vector that stores the kinematics chain
	// shared by fk, fk_jacobian, fk_hessian
	// shared by forwardkinematics.cpp and Digit_model_floating.cpp
	std::vector<int> chain;  

	// internal copies
	int current_order = -1;
	int end_copy = -1;
	int start_copy = -1;
	VecX q_copy;
	Transform T_copy;
	std::vector<Transform> dTdq_copy;
	std::vector<std::vector<Transform>> ddTddq_copy;
	std::vector<std::vector<std::vector<Transform>>> dddTdddq_copy;
};

}  // namespace IDTO

#endif