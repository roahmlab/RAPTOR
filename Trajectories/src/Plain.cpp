#include "Plain.h"

namespace IDTO {

Plain::Plain(const int Nact_input) {
    N = 1;
    Nact = Nact_input;
    varLength = Nact_input;

    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);
}

void Plain::compute(const VecX& z, bool compute_derivatives) {
    if (z.size() < varLength) {
        throw std::invalid_argument("Plain: decision variable vector has wrong size");
    }

    for (int i = 0; i < N; i++) {
        q(i) = z;
        // q_d(i) = 0;
        // q_dd(i) = 0;
    }

    if (compute_derivatives) {
        for (int i = 0; i < N; i++) {
            pq_pz(i).resize(Nact, varLength);
            pq_pz(i).setIdentity();
            // pq_d_pz(i).resize(Nact, varLength);
            // pq_dd_pz(i).resize(Nact, varLength);
        }
    }
}

}; // namespace IDTO

