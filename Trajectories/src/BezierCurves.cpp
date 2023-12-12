#include "BezierCurves.h"

namespace IDTO {

// UN FINISHED

BezierCurves::BezierCurves(const VecX& tspan_input, int Nact_input, int degree_input) : 
    Trajectories(tspan_input, Nact_input),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);
}

BezierCurves::BezierCurves(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input) :
    Trajectories(T_input, N_input, Nact_input, time_discretization),
    degree(degree_input) {
    varLength = (degree + 1) * Nact;

    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);
}

void BezierCurves::compute(const VecX& z, bool compute_derivatives) {
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

        B(0) = 1;
        dB(0) = 0;
        ddB(0) = 0;

        B(1) = t;
        dB(1) = 1;
        ddB(1) = 0;

        for (int i = 2; i <= degree; i++) {
            B(i) = pow(t, i) / (i * (i-1));
            dB(i) = pow(t, i-1) / (i-1);
            ddB(i) = pow(t, i-2);
        }

        q(x) = coefficients.transpose() * B;
        q_d(x) = coefficients.transpose() * dB;
        q_dd(x) = coefficients.transpose() * ddB;

        if (compute_derivatives) {
            for (int i = 0; i < Nact; i++) {
                pq_pz(x).block(i, i * (degree + 1), 1, degree + 1) = B.transpose();
                pq_d_pz(x).block(i, i * (degree + 1), 1, degree + 1) = dB.transpose();
                pq_dd_pz(x).block(i, i * (degree + 1), 1, degree + 1) = ddB.transpose();
            }
        }
    }
}

}; // namespace IDTO

