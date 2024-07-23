#include "Plain.h"

namespace RAPTOR {

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
}

Plain::Plain(const int N_input, 
             const int Nact_input) {
    N = N_input;
    Nact = Nact_input;
    varLength = Nact_input * N_input;

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
}

void Plain::compute(const VecX& z, 
                    bool compute_derivatives,
                    bool compute_hessian) {
    if (z.size() < varLength) {
        std::cerr << "function input: z.size() = " << z.size() << std::endl;
        std::cerr << "desired: varLength = " << varLength << std::endl;
        throw std::invalid_argument("Plain: decision variable vector has wrong size");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    for (int i = 0; i < N; i++) {
        q(i) = z.block(i * Nact, 0, Nact, 1);
        // q_d(i) = 0;
        // q_dd(i) = 0;

        if (compute_derivatives) {
            pq_pz(i).resize(Nact, varLength);
            pq_pz(i).block(0, i * Nact, Nact, Nact) = Eigen::MatrixXd::Identity(Nact, Nact);
            // pq_d_pz(i).resize(Nact, varLength);
            // pq_dd_pz(i).resize(Nact, varLength);
        }

        if (compute_hessian) {
            for (int j = 0; j < Nact; j++) {
                pq_pz_pz(j, i).setZero();
                // pq_d_pz_pz(j, i).setZero();
                // pq_dd_pz_pz(j, i).setZero();
            }
        }
    }
}

}; // namespace RAPTOR

