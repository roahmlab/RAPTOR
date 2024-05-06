#include "TransferBezierCurves.h"

namespace IDTO {

TransferBezierCurves::TransferBezierCurves(double T_input, 
                                           int N_input, 
                                           int Nact_input, 
                                           const int degree_input,
                                           TimeDiscretization time_discretization, 
                                           const TransferTrajectoryParameters& ttp_input) {
    T = T_input;
    N = N_input;
    Nact = Nact_input;
    degree = degree_input;
    ttp = ttp_input;
    
    tspan = VecX::Zero(N);
    for (int i = 0; i < N; i++) {
        tspan(i) = 0.5 * T * (1 - cos(M_PI * (2 * i + 1) / (2 * N)));
    }

    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);

    Bionomials = VecX::Ones(degree + 1);
    for (int j = 1; j <= degree / 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (degree + 1 - j) / j;
        Bionomials(degree - j) = Bionomials(j);
    }
    
    if (ttp.q0.size() != Nact) {
        throw std::invalid_argument("TransferBezierCurves: q0.size() != Nact");
    }
    if (ttp.qT.size() != Nact) {
        throw std::invalid_argument("TransferBezierCurves: qT.size() != Nact");
    }

    varLength = Nact * (degree - 5);

    coefficients.resize(degree + 1, Nact);
    coefficients.row(0) = ttp.q0.transpose();
    coefficients.row(1) = coefficients.row(0);
    coefficients.row(2) = coefficients.row(0);
    coefficients.row(degree) = ttp.qT.transpose();
    coefficients.row(degree - 1) = coefficients.row(degree);
    coefficients.row(degree - 2) = coefficients.row(degree);

    if (varLength <= 0) {
        throw std::invalid_argument("TransferBezierCurves: degree must be larget than 5");
    }

    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);
}

void TransferBezierCurves::compute(const VecX& z, bool compute_derivatives) {
    if (z.size() != varLength) {
        throw std::invalid_argument("TransferBezierCurves: decision variable vector has wrong size");
    }

    if (if_computed(z, compute_derivatives)) return;

    for (int i = 3; i < degree - 2; i++) {
        coefficients.row(i) = z.block((i - 3) * Nact, 0, Nact, 1).transpose();
    }

    for (int x = 0; x < N; x++) {
        double t = tspan(x) / T;

        q(x) = VecX::Zero(Nact + 1);
        q_d(x) = VecX::Zero(Nact + 1);
        q_dd(x) = VecX::Zero(Nact + 1);

        if (compute_derivatives) {
            pq_pz(x) = MatX::Zero(Nact, varLength);
            pq_d_pz(x) = MatX::Zero(Nact, varLength);
            pq_dd_pz(x) = MatX::Zero(Nact, varLength);
        }

        // Compute tA(i, j) = t(i)^j, 
        //         tB(i, j) = (1-t(i))^(degree-j)
        VecX tA = VecX::Ones(degree + 1);
        VecX tB = VecX::Ones(degree + 1);
        VecX dtA = VecX::Zero(degree + 1);
        VecX dtB = VecX::Zero(degree + 1);
        VecX ddtA = VecX::Zero(degree + 1);
        VecX ddtB = VecX::Zero(degree + 1);

        // Loop to compute tA and tB
        for (int j = 1; j <= degree; j++) {
            tA(j) = t * tA(j - 1);
            tB(degree - j) = (1 - t) * tB(degree - j + 1);

            dtA(j) = j * tA(j - 1);
            dtB(degree - j) = -j * tB(degree - j + 1);

            ddtA(j) = j * dtA(j - 1);
            ddtB(degree - j) = -j * dtB(degree - j + 1);
        }

        B = Bionomials.array() * tA.array() * tB.array();
        dB = Bionomials.array() * (dtA.array() * tB.array() +  
                                   tA.array() * dtB.array()) / T;
        ddB = Bionomials.array() * (ddtA.array() * tB.array() + 
                                    2 * dtA.array() * dtB.array() + 
                                    tA.array() * ddtB.array()) / (T * T);

        q(x).head(Nact)    = coefficients.transpose() * B;
        q_d(x).head(Nact)  = coefficients.transpose() * dB;
        q_dd(x).head(Nact) = coefficients.transpose() * ddB;

        if (compute_derivatives) {
            for (int i = 0; i < Nact; i++) {
                for (int j = 3; j < degree - 2; j++) {
                    pq_pz(x)(i, (j - 3) * Nact + i) = B(j);
                    pq_d_pz(x)(i, (j - 3) * Nact + i) = dB(j);
                    pq_dd_pz(x)(i, (j - 3) * Nact + i) = ddB(j);
                }
            }
        }
    }
}

}; // namespace IDTO