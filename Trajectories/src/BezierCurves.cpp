#include "BezierCurves.h"

namespace IDTO {

BezierCurves::BezierCurves(const VecX& tspan_input, int Nact_input, int degree_input) : 
    Trajectories(tspan_input, Nact_input),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);

    Bionomials = VecX::Ones(degree + 1);
    for (int j = 1; j <= degree / 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (degree + 1 - j) / j;
        Bionomials(degree - j) = Bionomials(j);
    }
}

BezierCurves::BezierCurves(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input) :
    Trajectories(T_input, N_input, Nact_input, time_discretization),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);

    Bionomials = VecX::Ones(degree + 1);
    for (int j = 1; j <= degree / 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (degree + 1 - j) / j;
        Bionomials(degree - j) = Bionomials(j);
    }
}

void BezierCurves::compute(const VecX& z, bool compute_derivatives) {
    if (z.size() < varLength) {
        throw std::invalid_argument("BezierCurves: decision variable vector has wrong size");
    }

    MatX coefficients = z.head((degree + 1) * Nact).reshaped(degree + 1, Nact);

    for (int x = 0; x < N; x++) {
        double t = tspan(x) / T;

        q(x) = VecX::Zero(Nact);
        q_d(x) = VecX::Zero(Nact);
        q_dd(x) = VecX::Zero(Nact);

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

        q(x)    = coefficients.transpose() * B;
        q_d(x)  = coefficients.transpose() * dB;
        q_dd(x) = coefficients.transpose() * ddB;

        if (compute_derivatives) {
            for (int i = 0; i < Nact; i++) {
                pq_pz(x).block(i, i * (degree + 1), 1, degree + 1)    = B.transpose();
                pq_d_pz(x).block(i, i * (degree + 1), 1, degree + 1)  = dB.transpose();
                pq_dd_pz(x).block(i, i * (degree + 1), 1, degree + 1) = ddB.transpose();
            }
        }
    }
}

}; // namespace IDTO

