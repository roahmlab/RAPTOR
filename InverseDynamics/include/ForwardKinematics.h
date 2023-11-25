#ifndef FORWARD_KINEMATICS_HPP
#define FORWARD_KINEMATICS_HPP

#include "Transform.h"

namespace IDTO {

class ForwardKinematicsHighOrderDerivative {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;

    // Constructor
	ForwardKinematicsHighOrderDerivative() = default;

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

	VecX Transform2xyzrpy(const Transform& T);

	void Transform2xyzrpyJacobian(Eigen::MatrixXd& J, 
								  const Transform& T, 
								  const std::vector<Transform>& dTdq);

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
};

}  // namespace IDTO

#endif