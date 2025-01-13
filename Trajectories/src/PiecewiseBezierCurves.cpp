#include "PiecewiseBezierCurves.h"

namespace RAPTOR {

PiecewiseBezierCurves::PiecewiseBezierCurves(double T_input, 
                                             int N_input, 
                                             int Nact_input, 
                                             int degree_input,
                                             const VecX q0_input,
                                             const VecX qT_input) {
    T = T_input;
    N = N_input;
    Nact = Nact_input;
    degree = degree_input;
    q0 = q0_input;
    qT = qT_input;

    if (N % (degree + 1) != 0) {
        throw std::invalid_argument("PiecewiseBezierCurves: N must be a multiple of (degree + 1)");
    }
    
    tspan = VecX::LinSpaced(N, 0, T);

    B = VecX::Zero(6);
    dB = VecX::Zero(6);
    ddB = VecX::Zero(6);

    Bionomials = VecX::Ones(6);
    for (int j = 1; j <= 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (6 - j) / j;
        Bionomials(5 - j) = Bionomials(j);
    }

    varLength = Nact * degree * 3;
    
    if (q0.size() != Nact) {
        optimize_begin_position = true;
        begin_position_offset = varLength;
        varLength += Nact;
    }
    if (qT.size() != Nact) {
        optimize_end_position = true;
        end_position_offset = varLength;
        varLength += Nact;
    }

    coefficients.resize(6, Nact);

    initialize_memory();
}

void PiecewiseBezierCurves::compute(const VecX& z, 
                                    bool compute_derivatives,
                                    bool compute_hessian) {
    if (z.size() != varLength) {
        std::cerr << "function input: z.size() = " << z.size() << std::endl;
        std::cerr << "desired: varLength = " << varLength << std::endl;
        throw std::invalid_argument("PiecewiseBezierCurves: decision variable vector has wrong size");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    const int seg_size = N / (degree + 1);
    const double localT = T / (degree + 1);

    for (int x = 0; x < N; x++) {
        // determine which Bezier curve we are at right now
        int id = x / seg_size;

        // determine the local time within the Bezier curve (scaled to [0, 1])
        double t = (tspan(x) - id * localT) / localT;

        // useful constants
        double velocity_factor = localT / 5.0;
        double acceleration_factor = localT * localT / 20.0;

        // determine the initial and end position, velocity and acceleration for each piece of Bezier curve
        VecX begin_q, begin_q_d, begin_q_dd;
        VecX end_q, end_q_d, end_q_dd;

        if (id == 0) {
            const VecX& control_point = z.head(3 * Nact);

            if (optimize_begin_position) {
                begin_q = z.segment(begin_position_offset, Nact).transpose();
            }
            else {
                begin_q = q0.transpose();
            }

            begin_q_d = VecX::Zero(Nact);
            begin_q_dd = VecX::Zero(Nact);
            
            end_q = control_point.head(Nact).transpose();
            end_q_d = control_point.segment(Nact, Nact).transpose();
            end_q_dd = control_point.tail(Nact).transpose();
        }
        else if (id == degree) {
            const VecX& control_point = z.segment((id - 1) * 3 * Nact, 3 * Nact);

            begin_q = control_point.head(Nact).transpose();
            begin_q_d = control_point.segment(Nact, Nact).transpose();
            begin_q_dd = control_point.tail(Nact).transpose();

            if (optimize_end_position) {
                end_q = z.segment(end_position_offset, Nact).transpose();
            }
            else {
                end_q = qT.transpose();
            }
            end_q_d = VecX::Zero(Nact);
            end_q_dd = VecX::Zero(Nact);
        }
        else {
            const VecX& control_point1 = z.segment((id - 1) * 3 * Nact, 3 * Nact);
            const VecX& control_point2 = z.segment(id * 3 * Nact, 3 * Nact);
            begin_q = control_point1.head(Nact).transpose();
            begin_q_d = control_point1.segment(Nact, Nact).transpose();
            begin_q_dd = control_point1.tail(Nact).transpose();
            end_q = control_point2.head(Nact).transpose();
            end_q_d = control_point2.segment(Nact, Nact).transpose();
            end_q_dd = control_point2.tail(Nact).transpose();
        }   

        // compute the Bezier coefficients 
        // based on the initial and end position, velocity and acceleration
        coefficients.row(0) = begin_q;
        coefficients.row(1) = begin_q + 
                              begin_q_d * velocity_factor;
        coefficients.row(2) = begin_q + 
                              begin_q_d * 2.0 * velocity_factor + 
                              begin_q_dd * acceleration_factor;
        coefficients.row(3) = end_q - 
                              end_q_d * 2.0 * velocity_factor + 
                              end_q_dd * acceleration_factor;
        coefficients.row(4) = end_q - 
                              end_q_d * velocity_factor;
        coefficients.row(5) = end_q;

        // now compute the time-related terms (Bernstein polynomial basis) in the Bezier curve
        // Compute tA(i, j) = t(i)^j, 
        //         tB(i, j) = (1-t(i))^(degree-j)
        VecX tA = VecX::Ones(6);
        VecX tB = VecX::Ones(6);
        VecX dtA = VecX::Zero(6);
        VecX dtB = VecX::Zero(6);
        VecX ddtA = VecX::Zero(6);
        VecX ddtB = VecX::Zero(6);

        // Loop to compute tA and tB
        for (int j = 1; j <= 5; j++) {
            tA(j) = t * tA(j - 1);
            tB(5 - j) = (1 - t) * tB(5 - j + 1);

            dtA(j) = j * tA(j - 1);
            dtB(5 - j) = -j * tB(5 - j + 1);

            ddtA(j) = j * dtA(j - 1);
            ddtB(5 - j) = -j * dtB(5 - j + 1);
        }

        B = Bionomials.array() * tA.array() * tB.array();
        dB = Bionomials.array() * (dtA.array() * tB.array() +  
                                   tA.array() * dtB.array()) / localT;
        ddB = Bionomials.array() * (ddtA.array() * tB.array() + 
                                    2 * dtA.array() * dtB.array() + 
                                    tA.array() * ddtB.array()) / (localT * localT);

        // finally, compute the trajectory by multiplying the Bezier coefficients 
        // with the time-related terms (Bernstein polynomial basis)
        q(x) = VecX::Zero(Nact);
        q_d(x) = VecX::Zero(Nact);
        q_dd(x) = VecX::Zero(Nact);

        q(x)    = coefficients.transpose() * B;
        q_d(x)  = coefficients.transpose() * dB;
        q_dd(x) = coefficients.transpose() * ddB;

        if (compute_derivatives) {
            pq_pz(x) = MatX::Zero(Nact, varLength);
            pq_d_pz(x) = MatX::Zero(Nact, varLength);
            pq_dd_pz(x) = MatX::Zero(Nact, varLength);

            if (id == 0) {
                for (int i = 0; i < Nact; i++) {
                    pq_pz(x)(i, i) = B(3) + B(4) + B(5);
                    pq_pz(x)(i, i + Nact) = -B(3) * 2.0 * velocity_factor - B(4) * velocity_factor;
                    pq_pz(x)(i, i + 2 * Nact) = B(3) * acceleration_factor;
                    pq_d_pz(x)(i, i) = dB(3) + dB(4) + dB(5);
                    pq_d_pz(x)(i, i + Nact) = -dB(3) * 2.0 * velocity_factor - dB(4) * velocity_factor;
                    pq_d_pz(x)(i, i + 2 * Nact) = dB(3) * acceleration_factor;
                    pq_dd_pz(x)(i, i) = ddB(3) + ddB(4) + ddB(5);
                    pq_dd_pz(x)(i, i + Nact) = -ddB(3) * 2.0 * velocity_factor - ddB(4) * velocity_factor;
                    pq_dd_pz(x)(i, i + 2 * Nact) = ddB(3) * acceleration_factor;

                    if (optimize_begin_position) {
                        pq_pz(x)(i, begin_position_offset + i) = B(0) + B(1) + B(2);
                        pq_d_pz(x)(i, begin_position_offset + i) = dB(0) + dB(1) + dB(2);
                        pq_dd_pz(x)(i, begin_position_offset + i) = ddB(0) + ddB(1) + ddB(2);
                    }
                }
            }
            else if (id == degree) {
                for (int i = 0; i < Nact; i++) {
                    pq_pz(x)(i, (degree - 1) * 3 * Nact + i) = B(0) + B(1) + B(2);
                    pq_pz(x)(i, (degree - 1) * 3 * Nact + i + Nact) = B(1) * velocity_factor + B(2) * 2.0 * velocity_factor;
                    pq_pz(x)(i, (degree - 1) * 3 * Nact + i + 2 * Nact) = B(2) * acceleration_factor;
                    pq_d_pz(x)(i, (degree - 1) * 3 * Nact + i) = dB(0) + dB(1) + dB(2);
                    pq_d_pz(x)(i, (degree - 1) * 3 * Nact + i + Nact) = dB(1) * velocity_factor + dB(2) * 2.0 * velocity_factor;
                    pq_d_pz(x)(i, (degree - 1) * 3 * Nact + i + 2 * Nact) = dB(2) * acceleration_factor;
                    pq_dd_pz(x)(i, (degree - 1) * 3 * Nact + i) = ddB(0) + ddB(1) + ddB(2);
                    pq_dd_pz(x)(i, (degree - 1) * 3 * Nact + i + Nact) = ddB(1) * velocity_factor + ddB(2) * 2.0 * velocity_factor;
                    pq_dd_pz(x)(i, (degree - 1) * 3 * Nact + i + 2 * Nact) = ddB(2) * acceleration_factor;

                    if (optimize_end_position) {
                        pq_pz(x)(i, end_position_offset + i) = B(3) + B(4) + B(5);
                        pq_d_pz(x)(i, end_position_offset + i) = dB(3) + dB(4) + dB(5);
                        pq_dd_pz(x)(i, end_position_offset + i) = ddB(3) + ddB(4) + ddB(5);
                    }
                }
            }
            else {
                for (int i = 0; i < Nact; i++) {
                    pq_pz(x)(i, (id - 1) * 3 * Nact + i) = B(0) + B(1) + B(2);
                    pq_pz(x)(i, (id - 1) * 3 * Nact + i + Nact) = B(1) * velocity_factor + B(2) * 2.0 * velocity_factor;
                    pq_pz(x)(i, (id - 1) * 3 * Nact + i + 2 * Nact) = B(2) * acceleration_factor;
                    pq_d_pz(x)(i, (id - 1) * 3 * Nact + i) = dB(0) + dB(1) + dB(2);
                    pq_d_pz(x)(i, (id - 1) * 3 * Nact + i + Nact) = dB(1) * velocity_factor + dB(2) * 2.0 * velocity_factor;
                    pq_d_pz(x)(i, (id - 1) * 3 * Nact + i + 2 * Nact) = dB(2) * acceleration_factor;
                    pq_dd_pz(x)(i, (id - 1) * 3 * Nact + i) = ddB(0) + ddB(1) + ddB(2);
                    pq_dd_pz(x)(i, (id - 1) * 3 * Nact + i + Nact) = ddB(1) * velocity_factor + ddB(2) * 2.0 * velocity_factor;
                    pq_dd_pz(x)(i, (id - 1) * 3 * Nact + i + 2 * Nact) = ddB(2) * acceleration_factor;

                    pq_pz(x)(i, id * 3 * Nact + i) = B(3) + B(4) + B(5);
                    pq_pz(x)(i, id * 3 * Nact + i + Nact) = -B(3) * 2.0 * velocity_factor - B(4) * velocity_factor;
                    pq_pz(x)(i, id * 3 * Nact + i + 2 * Nact) = B(3) * acceleration_factor;
                    pq_d_pz(x)(i, id * 3 * Nact + i) = dB(3) + dB(4) + dB(5);
                    pq_d_pz(x)(i, id * 3 * Nact + i + Nact) = -dB(3) * 2.0 * velocity_factor - dB(4) * velocity_factor;
                    pq_d_pz(x)(i, id * 3 * Nact + i + 2 * Nact) = dB(3) * acceleration_factor;
                    pq_dd_pz(x)(i, id * 3 * Nact + i) = ddB(3) + ddB(4) + ddB(5);
                    pq_dd_pz(x)(i, id * 3 * Nact + i + Nact) = -ddB(3) * 2.0 * velocity_factor - ddB(4) * velocity_factor;
                    pq_dd_pz(x)(i, id * 3 * Nact + i + 2 * Nact) = ddB(3) * acceleration_factor;
                }
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