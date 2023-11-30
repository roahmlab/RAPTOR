#include "Polynomials.h"

namespace IDTO {

Polynomials::Polynomials(const VecX& tspan_input, int Nact_input, int degree_input) : 
    Trajectories(tspan_input, Nact_input),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    P = VecX::Zero(degree + 1);
    dP = VecX::Zero(degree + 1);
    ddP = VecX::Zero(degree + 1);
}

Polynomials::Polynomials(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input) :
    Trajectories(T_input, N_input, Nact_input, time_discretization),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    P = VecX::Zero(degree + 1);
    dP = VecX::Zero(degree + 1);
    ddP = VecX::Zero(degree + 1);
}

void Polynomials::compute(const VecX& z, bool compute_derivatives) {
    assert(z.size() == varLength);

    MatX coefficients = z.reshaped(degree + 1, Nact);

    for (int x = 0; x < N; x++) {
        double t = tspan(x);

        q(x) = VecX::Zero(Nact);
        q_d(x) = VecX::Zero(Nact);
        q_dd(x) = VecX::Zero(Nact);

        if (compute_derivatives) {
            pq_pz(x) = MatX::Zero(Nact, varLength);
            pq_d_pz(x) = MatX::Zero(Nact, varLength);
            pq_dd_pz(x) = MatX::Zero(Nact, varLength);
        }

        P(0) = 1;
        dP(0) = 0;
        ddP(0) = 0;

        P(1) = t;
        dP(1) = 1;
        ddP(1) = 0;

        for (int i = 2; i <= degree; i++) {
            P(i) = pow(t, i) / (i * (i-1));
            dP(i) = pow(t, i-1) / (i-1);
            ddP(i) = pow(t, i-2);
        }

        q(x) = coefficients.transpose() * P;
        q_d(x) = coefficients.transpose() * dP;
        q_dd(x) = coefficients.transpose() * ddP;

        if (compute_derivatives) {
            for (int i = 0; i < Nact; i++) {
                pq_pz(x).block(i, i * (degree + 1), 1, degree + 1) = P.transpose();
                pq_d_pz(x).block(i, i * (degree + 1), 1, degree + 1) = dP.transpose();
                pq_dd_pz(x).block(i, i * (degree + 1), 1, degree + 1) = ddP.transpose();
            }
        }
    }
}

}; // namespace IDTO

