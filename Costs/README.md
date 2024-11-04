# Costs

This class implements different types of cost functions that will be later put into the optimization.
Specifically, we would like to provide the formulation `f(z)`.
The core function is called `compute(const Eigen::VectorXd& z, bool compute_derivatives, bool compute_hessian`. 
The cost class usually requires a shared pointer of a trajectory class in the constructor, or a shared pointer of a forward kinematics and inverse dynamics class depending on the purpose of the cost functions. 
The results are stored in `f`, which is a class member, for the optimizer class to access.
The gradient is stored in `grad_f` and the hessian is stored in `hess_f`

