#include "WaitrBezierCurves.h"

namespace RAPTOR {

WaitrBezierCurves::WaitrBezierCurves(const VecX& tspan_input, 
                                     int Nact_input, 
                                     const WaitrTrajectoryParameters& atp_input) : 
    BezierCurves(tspan_input, Nact_input, WAITR_BEZIER_CURVE_DEGREE),
    atp(atp_input) {
    if (atp.q0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q0.size() != Nact");
    }
    if (atp.q_d0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q_d0.size() != Nact");
    }
    if (atp.q_dd0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q_dd0.size() != Nact");
    }

    varLength = Nact;

    coefficients.resize(degree + 1, Nact);
    coefficients.row(0) = atp.q0.transpose();
    coefficients.row(1) = coefficients.row(0) + 
                          (T * atp.q_d0.transpose()) / 5.0;
    coefficients.row(2) = coefficients.row(0) + 
                          (2.0 * T * atp.q_d0.transpose()) / 5.0 + 
                          (T * T * atp.q_dd0.transpose()) / 20.0;
}

WaitrBezierCurves::WaitrBezierCurves(double T_input, 
                                     int N_input, 
                                     int Nact_input, 
                                     TimeDiscretization time_discretization, 
                                     const WaitrTrajectoryParameters& atp_input) :
    BezierCurves(T_input, N_input, Nact_input, time_discretization, WAITR_BEZIER_CURVE_DEGREE),
    atp(atp_input) {
    if (atp.q0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q0.size() != Nact");
    }
    if (atp.q_d0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q_d0.size() != Nact");
    }
    if (atp.q_dd0.size() != Nact) {
        throw std::invalid_argument("WaitrBezierCurves: q_dd0.size() != Nact");
    }

    varLength = Nact;

    coefficients.resize(degree + 1, Nact);
    coefficients.row(0) = atp.q0.transpose();
    coefficients.row(1) = atp.q0.transpose() + 
                          (T * atp.q_d0.transpose()) / 5.0;
    coefficients.row(2) = atp.q0.transpose() + 
                          (2.0 * T * atp.q_d0.transpose()) / 5.0 + 
                          (T * T * atp.q_dd0.transpose()) / 20.0;
}

void WaitrBezierCurves::compute(const VecX& z, 
                                bool compute_derivatives,
                                bool compute_hessian) {
    if (z.size() != varLength) {
        std::cerr << "function input: z.size() = " << z.size() << std::endl;
        std::cerr << "desired: varLength = " << varLength << std::endl;
        throw std::invalid_argument("WaitrBezierCurves: decision variable vector has wrong size");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    coefficients.row(5) = z.transpose();
    coefficients.row(4) = z.transpose();
    coefficients.row(3) = z.transpose();

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
                pq_pz(x)(i, i)    = B(3) + B(4) + B(5);
                pq_d_pz(x)(i, i)  = dB(3) + dB(4) + dB(5);
                pq_dd_pz(x)(i, i) = ddB(3) + ddB(4) + ddB(5);
            }
        }

        if (compute_hessian) {
            for (int i = 0; i < Nact; i++) {
                pq_pz_pz(i, x).setZero();
                pq_d_pz_pz(i, x).setZero();
                pq_dd_pz_pz(i, x).setZero();
            }
        }
    }
}

}; // namespace RAPTOR