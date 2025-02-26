#include "BezierCurves.h"

namespace RAPTOR {

BezierCurves::BezierCurves(const VecX& tspan_input, 
                           int Nact_input, 
                           int degree_input) : 
    Trajectories((degree_input + 1) * Nact_input, tspan_input, Nact_input),
    degree(degree_input) {
    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);

    Bionomials = VecX::Ones(degree + 1);
    for (int j = 1; j <= degree / 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (degree + 1 - j) / j;
        Bionomials(degree - j) = Bionomials(j);
    }
}

BezierCurves::BezierCurves(double T_input, 
                           int N_input, 
                           int Nact_input, 
                           TimeDiscretization time_discretization, 
                           int degree_input) :
    Trajectories((degree_input + 1) * Nact_input, T_input, N_input, Nact_input, time_discretization),
    degree(degree_input) {
    B = VecX::Zero(degree + 1);
    dB = VecX::Zero(degree + 1);
    ddB = VecX::Zero(degree + 1);

    Bionomials = VecX::Ones(degree + 1);
    for (int j = 1; j <= degree / 2; j++) {
        Bionomials(j) = Bionomials(j - 1) * (degree + 1 - j) / j;
        Bionomials(degree - j) = Bionomials(j);
    }
}

void BezierCurves::constrainInitialPosition(const VecX& q0_input) {
    if (q0_input.size() != Nact) {
        std::cerr << "function input: q0_input.size() = " << q0_input.size() << std::endl;
        std::cerr << "desired: Nact = " << Nact << std::endl;
        throw std::invalid_argument("BezierCurves: initial position vector has wrong size");
    }

    q0 = q0_input;
    constrain_initial_position = true;

    varLength -= Nact;
    initialize_memory();
}

void BezierCurves::constrainInitialVelocity(const VecX& q_d0_input) {
    if (!constrain_initial_position) {
        throw std::invalid_argument("BezierCurves: initial position must be constrained first");
    }

    if (q_d0_input.size() != Nact) {
        std::cerr << "function input: q_d0_input.size() = " << q_d0_input.size() << std::endl;
        std::cerr << "desired: Nact = " << Nact << std::endl;
        throw std::invalid_argument("BezierCurves: initial velocity vector has wrong size");
    }

    q_d0 = q_d0_input;
    constrain_initial_velocity = true;

    varLength -= Nact;
    initialize_memory();
}

void BezierCurves::compute(const VecX& z, 
                           bool compute_derivatives,
                           bool compute_hessian) {
    if (z.size() < varLength) {
        std::cerr << "function input: z.size() = " << z.size() << std::endl;
        std::cerr << "desired: varLength = " << varLength << std::endl;
        throw std::invalid_argument("BezierCurves: decision variable vector has wrong size");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    MatX coefficients = MatX::Zero(degree + 1, Nact);
    
    if (constrain_initial_position) {
        if (constrain_initial_velocity) {
            coefficients.topRows(1) = q0.transpose();
            coefficients.middleRows(1, 1) = q0.transpose() + q_d0.transpose() * T / 5.0;
            coefficients.bottomRows(degree - 1) = Utils::reshape(z.head((degree - 1) * Nact), 
                                                                 degree - 1, Nact);
        }
        else {
            coefficients.topRows(1) = q0.transpose();
            coefficients.bottomRows(degree) = Utils::reshape(z.head(degree * Nact), 
                                                             degree, Nact);
        }
    }
    else {
        coefficients = Utils::reshape(z.head((degree + 1) * Nact), 
                                      degree + 1, Nact);
    }

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
            if (constrain_initial_position) {
                if (constrain_initial_velocity) {
                    for (int i = 0; i < Nact; i++) {
                        pq_pz(x).block(i, i * (degree - 1), 1, degree - 1)    = B.tail(degree - 1).transpose();
                        pq_d_pz(x).block(i, i * (degree - 1), 1, degree - 1)  = dB.tail(degree - 1).transpose();
                        pq_dd_pz(x).block(i, i * (degree - 1), 1, degree - 1) = ddB.tail(degree - 1).transpose();
                    }
                }
                else {
                    for (int i = 0; i < Nact; i++) {
                        pq_pz(x).block(i, i * degree, 1, degree)    = B.tail(degree).transpose();
                        pq_d_pz(x).block(i, i * degree, 1, degree)  = dB.tail(degree).transpose();
                        pq_dd_pz(x).block(i, i * degree, 1, degree) = ddB.tail(degree).transpose();
                    }
                }
            }
            else {
                for (int i = 0; i < Nact; i++) {
                    pq_pz(x).block(i, i * (degree + 1), 1, degree + 1)    = B.transpose();
                    pq_d_pz(x).block(i, i * (degree + 1), 1, degree + 1)  = dB.transpose();
                    pq_dd_pz(x).block(i, i * (degree + 1), 1, degree + 1) = ddB.transpose();
                }
            }
        }

        if (compute_hessian) {
            // this has been done in the constructor
            // for (int i = 0; i < Nact; i++) {
            //     pq_pz_pz(i, x).setZero(); 
            //     pq_d_pz_pz(i, x).setZero(); 
            //     pq_dd_pz_pz(i, x).setZero(); 
            // }
        }
    }
}

}; // namespace RAPTOR

