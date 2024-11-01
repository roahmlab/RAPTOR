#include "Costs.h"

namespace RAPTOR {

Costs::Costs(int varLength) {
    initialize_memory(varLength);
}

bool Costs::is_computed(const VecX& z, 
                        bool compute_derivatives,
                        bool compute_hessian) {
    if (compute_hessian && !compute_derivatives) {
        throw std::invalid_argument("compute_derivatives needs to be true when compute_hessian is true.");
        return false;
    }

    if (!Utils::ifTwoVectorEqual(current_z, z, 0)) {
        current_z = z;
        if_compute_derivatives = compute_derivatives;
        if_compute_hessian = compute_hessian;
        return false;
    }

    if (compute_derivatives != if_compute_derivatives) {
        current_z = z;
        if_compute_derivatives = compute_derivatives;
        if_compute_hessian = compute_hessian;
        return false;
    }

    // current_z = z;  
    if_compute_derivatives = compute_derivatives;
    if_compute_hessian = compute_hessian;
    return true;
}

void Costs::initialize_memory(const int varLength) {
    f = 0;
    grad_f = VecX::Zero(varLength);
    hess_f = MatX::Zero(varLength, varLength);
}

}; // namespace RAPTOR