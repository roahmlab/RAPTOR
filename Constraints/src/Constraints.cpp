#include "Constraints.h"

namespace IDTO {

Constraints::Constraints(int m_input,
                         int varLength) {
    initialize_memory(m_input, varLength);
}

bool Constraints::is_computed(const VecX& z, 
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

void Constraints::initialize_memory(const int m_input, 
                                    const int varLength) {
    m = m_input;

    g = VecX::Zero(m);
    
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);

    pg_pz = MatX::Zero(m, varLength);

    pg_pz_pz.resize(m);
    for (int i = 0; i < m; i++) {
        pg_pz_pz(i) = MatX::Zero(varLength, varLength);
    }
}

}; // namespace IDTO