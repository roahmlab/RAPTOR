#include "Trajectories.h"

namespace IDTO {

Trajectories::Trajectories(const VecX& tspan_input, 
                           int Nact_input) :
    tspan(tspan_input),
    Nact(Nact_input){
    T = tspan.bottomLeftCorner<1,1>().value();
    N = tspan.size();

    varLength = 1;
    
    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);

    current_z.resize(1);
}

Trajectories::Trajectories(double T_input,  
                           int N_input, 
                           int Nact_input, 
                           TimeDiscretization time_discretization) :
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

    varLength = 1;
    
    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);
}

bool Trajectories::if_computed(const VecX& z, bool compute_derivatives) {
    if (!ifTwoVectorEqual(current_z, z, 0)) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        if (compute_derivatives != if_compute_derivatives) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        // current_z = z;  
        if_compute_derivatives = compute_derivatives;
        return true;
}

void Trajectories::compute(const VecX& z, bool compute_derivatives) {
    if (if_computed(z, compute_derivatives)) return;

    for (int i = 0; i < N; i++) {
        q(i) = VecX::Zero(Nact);
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);
    }

    if (compute_derivatives) {
        for (int i = 0; i < N; i++) {
            pq_pz(i).resize(Nact, varLength);
            pq_d_pz(i).resize(Nact, varLength);
            pq_dd_pz(i).resize(Nact, varLength);
        }
    }
}

}; // namespace IDTO

