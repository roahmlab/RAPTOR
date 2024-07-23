#include "Trajectories.h"

namespace RAPTOR {

Trajectories::Trajectories(const int varLength_input,
                           const VecX& tspan_input, 
                           int Nact_input) :
    varLength(varLength_input),
    tspan(tspan_input),
    Nact(Nact_input){
    T = tspan.bottomLeftCorner<1,1>().value();
    N = tspan.size();
    
    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);

    pq_pz_pz.resize(Nact, N);
    pq_d_pz_pz.resize(Nact, N);
    pq_dd_pz_pz.resize(Nact, N);

    for (int i = 0; i < N; i++) {
        q(i) = VecX::Zero(Nact);
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);

        pq_pz(i) = MatX::Zero(Nact, varLength);
        pq_d_pz(i) = MatX::Zero(Nact, varLength);
        pq_dd_pz(i) = MatX::Zero(Nact, varLength);

        for (int j = 0; j < Nact; j++) {
            pq_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_d_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_dd_pz_pz(j, i) = MatX::Zero(varLength, varLength);
        }
    }

    current_z.resize(1);
}

Trajectories::Trajectories(const int varLength_input,
                           double T_input, 
                           int N_input, 
                           int Nact_input, 
                           TimeDiscretization time_discretization) :
    varLength(varLength_input),
    T(T_input),
    N(N_input),
    Nact(Nact_input) {;
    if (time_discretization == Uniform) {
        tspan = VecX::LinSpaced(N, 0, T);
    } 
    else if (time_discretization == Chebyshev) {
        tspan = VecX::Zero(N);
        for (int i = 1; i < N - 1; i++) {
            tspan(i) = 0.5 * T * (1 - cos(M_PI * (2 * i - 1) / (2 * (N - 2))));
        }
        tspan(0) = 0;
        tspan(N - 1) = T;
    }
    
    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);

    pq_pz_pz.resize(Nact, N);
    pq_d_pz_pz.resize(Nact, N);
    pq_dd_pz_pz.resize(Nact, N);

    for (int i = 0; i < N; i++) {
        q(i) = VecX::Zero(Nact);
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);

        pq_pz(i) = MatX::Zero(Nact, varLength);
        pq_d_pz(i) = MatX::Zero(Nact, varLength);
        pq_dd_pz(i) = MatX::Zero(Nact, varLength);

        for (int j = 0; j < Nact; j++) {
            pq_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_d_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_dd_pz_pz(j, i) = MatX::Zero(varLength, varLength);
        }
    }

    current_z.resize(1);
}

bool Trajectories::is_computed(const VecX& z, 
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

    if (compute_hessian != if_compute_hessian) {
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

void Trajectories::compute(const VecX& z, 
                           bool compute_derivatives,
                           bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    for (int i = 0; i < N; i++) {
        q(i) = VecX::Zero(Nact);
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);
    }

    if (compute_derivatives) {
        for (int i = 0; i < N; i++) {
            pq_pz(i) = MatX::Zero(Nact, varLength);
            pq_d_pz(i) = MatX::Zero(Nact, varLength);
            pq_dd_pz(i) = MatX::Zero(Nact, varLength);
        }
    }

    if (compute_hessian) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < Nact; j++) {
                pq_pz_pz(j, i) = MatX::Zero(varLength, varLength);
                pq_d_pz_pz(j, i) = MatX::Zero(varLength, varLength);
                pq_dd_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            }
        }
    }
}

}; // namespace RAPTOR

