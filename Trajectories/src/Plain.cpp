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

    for (int i = 0; i < N; i++) {
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);

        pq_d_pz(i).resize(Nact, varLength);
        pq_dd_pz(i).resize(Nact, varLength);

        pq_d_pz(i).setZero();
        pq_dd_pz(i).setZero();
    }
}

Plain::Plain(const int N_input, const int Nact_input) {
    N = N_input;
    Nact = Nact_input;
    varLength = Nact_input * N_input;

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

    if (if_computed(z, compute_derivatives)) return;

    for (int i = 0; i < N; i++) {
        q(i) = z.block(i * Nact, 0, Nact, 1);
        // q_d(i) = 0;
        // q_dd(i) = 0;
    }

    if (compute_derivatives) {
        for (int i = 0; i < N; i++) {
            pq_pz(i).resize(Nact, varLength);
            pq_pz(i).block(0, i * Nact, Nact, Nact) = Eigen::MatrixXd::Identity(Nact, Nact);
            // pq_d_pz(i).resize(Nact, varLength);
            // pq_dd_pz(i).resize(Nact, varLength);
        }
    }
}

}; // namespace IDTO

