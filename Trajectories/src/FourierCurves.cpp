#include "FourierCurves.h"

namespace RAPTOR {

FourierCurves::FourierCurves(const VecX& tspan_input, 
                             int Nact_input, 
                             int degree_input,
                             VecX q0_input,
                             VecX q_d0_input) : 
    Trajectories((2 * degree_input + 2) * Nact_input, tspan_input, Nact_input),
    degree(degree_input),
    q0(q0_input),
    q_d0(q_d0_input) {
    F = VecX::Zero(2 * degree + 1);
    dF = VecX::Zero(2 * degree + 1);
    ddF = VecX::Zero(2 * degree + 1);

    F0 = VecX::Zero(2 * degree + 1);
    dF0 = VecX::Zero(2 * degree + 1);

    pF_pw = VecX::Zero(2 * degree + 1);
    pdF_pw = VecX::Zero(2 * degree + 1);
    pddF_pw = VecX::Zero(2 * degree + 1);

    pF0_pw = VecX::Zero(2 * degree + 1);
    pdF0_pw = VecX::Zero(2 * degree + 1);

    if (q0.size() == Nact) {
        optimize_initial_position = false;
    }
    else {
        optimize_initial_position = true;
        varLength += Nact;
    }
    
    if (q_d0.size() == Nact) {
        optimize_initial_velocity = false;
    }
    else {
        optimize_initial_velocity = true;
        varLength += Nact;
    }

    // varLength is changed so we have to reallocate the memory
    for (int i = 0; i < N; i++) {
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

FourierCurves::FourierCurves(double T_input, 
                             int N_input, 
                             int Nact_input, 
                             TimeDiscretization time_discretization,
                             int degree_input,
                             VecX q0_input,
                             VecX q_d0_input) :
    Trajectories((2 * degree_input + 4) * Nact_input, T_input, N_input, Nact_input, time_discretization),
    degree(degree_input),
    q0(q0_input),
    q_d0(q_d0_input) {
    F = VecX::Zero(2 * degree + 1);
    dF = VecX::Zero(2 * degree + 1);
    ddF = VecX::Zero(2 * degree + 1);

    F0 = VecX::Zero(2 * degree + 1);
    dF0 = VecX::Zero(2 * degree + 1);

    pF_pw = VecX::Zero(2 * degree + 1);
    pdF_pw = VecX::Zero(2 * degree + 1);
    pddF_pw = VecX::Zero(2 * degree + 1);

    pF0_pw = VecX::Zero(2 * degree + 1);
    pdF0_pw = VecX::Zero(2 * degree + 1);

    if (q0.size() == Nact) {
        optimize_initial_position = false;
    }
    else {
        optimize_initial_position = true;
        varLength += Nact;
    }
    
    if (q_d0.size() == Nact) {
        optimize_initial_velocity = false;
    }
    else {
        optimize_initial_velocity = true;
        varLength += Nact;
    }

    // varLength is changed so we have to reallocate the memory
    for (int i = 0; i < N; i++) {
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

void FourierCurves::compute(const VecX& z, 
                            bool compute_derivatives,
                            bool compute_hessian) {
    if (z.size() < varLength) {
        std::cerr << "function input: z.size() = " << z.size() << std::endl;
        std::cerr << "desired: varLength = " << varLength << std::endl;
        throw std::invalid_argument("FourierCurves: decision variable vector has wrong size");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) return;

    if (compute_hessian) {
        throw std::invalid_argument("FourierCurves: Hessian computation not implemented");
    }

    Eigen::MatrixXd temp = z.head((2 * degree + 2) * Nact);
    MatX coefficients = Utils::reshape(temp, 2 * degree + 2, Nact);

    if (optimize_initial_position) {
        q0 = z.block((2 * degree + 2) * Nact, 0, Nact, 1);
    }
    if (optimize_initial_velocity) {
        if (optimize_initial_position) {
            q_d0 = z.block((2 * degree + 2) * Nact + Nact, 0, Nact, 1);
        }
        else {
            q_d0 = z.block((2 * degree + 2) * Nact, 0, Nact, 1);
        }
    }

    for (int x = 0; x < N; x++) {
        double t = tspan(x);

        q(x) = VecX::Zero(Nact);
        q_d(x) = VecX::Zero(Nact);
        q_dd(x) = VecX::Zero(Nact);

        if (compute_derivatives) {
            pq_pz(x).resize(Nact, varLength);
            pq_d_pz(x).resize(Nact, varLength);
            pq_dd_pz(x).resize(Nact, varLength);
            pq_pz(x).setZero();
            pq_d_pz(x).setZero();
            pq_dd_pz(x).setZero();
        }

        for (int i = 0; i < Nact; i++) {
            const VecX& kernel = coefficients.block(0, i, 2 * degree + 1, 1);
            double w = coefficients(2 * degree + 1, i);

            ddF(0) = 1;
            dF(0)  = t;
            F(0)   = t * t * 0.5;

            dF0(0) = 0;
            F0(0) = 0;

            if (compute_derivatives) {
                pF_pw(0) = 0;
                pdF_pw(0) = 0;
                pddF_pw(0) = 0;

                pF0_pw(0) = 0;
                pdF0_pw(0) = 0;
            }

            for (int j = 0; j < degree; j++) {
                double jt = (j + 1) * t;
                double sinjwt = sin(w * jt);
                double cosjwt = cos(w * jt);

                ddF(2 * j + 1) = cosjwt;
                ddF(2 * j + 2) = sinjwt;

                double jw = (j + 1) * w;
                dF(2 * j + 1) = sinjwt / jw;
                dF(2 * j + 2) = -cosjwt / jw;
                dF0(2 * j + 2) = -1 / jw;

                double j2w2 = jw * jw;
                F(2 * j + 1) = -cosjwt / j2w2;
                F(2 * j + 2) = -sinjwt / j2w2;
                F0(2 * j + 1) = -1 / j2w2;

                if (compute_derivatives) {
                    pddF_pw(2 * j + 1) = -jt * sinjwt;
                    pddF_pw(2 * j + 2) = jt * cosjwt;

                    double jw2 = (j + 1) * w * w;
                    pdF_pw(2 * j + 1) = (t * cosjwt) / w - sinjwt / jw2;
                    pdF_pw(2 * j + 2) = cosjwt / jw2 + (t * sinjwt) / w;
                    pdF0_pw(2 * j + 2) = 1 / jw2;

                    double j2w3 = jw2 * (j + 1) * w;
                    pF_pw(2 * j + 1) = (2 * cosjwt + (j + 1) * w * t * sinjwt) / j2w3;
                    pF_pw(2 * j + 2) = (2 * sinjwt - (j + 1) * w * t * cosjwt) / j2w3;
                    pF0_pw(2 * j + 1) = 2 / j2w3;
                }
            }

            q_dd(x)(i) = ddF.dot(kernel);

            double q_d_raw = dF.dot(kernel);
            double q_d_raw0 = dF0.dot(kernel);
            q_d(x)(i) = q_d_raw + (q_d0(i) - q_d_raw0);

            double q_raw = F.dot(kernel) + (q_d0(i) - q_d_raw0) * t;
            double q_raw0 = F0.dot(kernel);
            q(x)(i) = q_raw + (q0(i) - q_raw0);

            if (compute_derivatives) {
                // pq_pz
                // derivative with respect to a_i
                for (int j = 0; j < 2 * degree + 1; j++) {
                    pq_pz(x)(i, i * (2 * degree + 2) + j) = F(j) - F0(j) - dF0(j) * tspan(x);
                }

                // derivative with respect to w
                pq_pz(x)(i, i * (2 * degree + 2) + (2 * degree + 1)) = kernel.dot(pF_pw - pF0_pw - pdF0_pw * tspan(x));

                // derivative with respect to q0
                if (optimize_initial_position) {
                    pq_pz(x)(i, Nact * (2 * degree + 2) + i) = 1; 
                }

                // derivative with respect to q_d0
                if (optimize_initial_velocity) {
                    if (optimize_initial_position) {
                        pq_pz(x)(i, Nact * (2 * degree + 2) + Nact + i) = tspan(x);   
                    }
                    else {
                        pq_pz(x)(i, Nact * (2 * degree + 2) + i) = tspan(x);
                    }
                }

                // pq_d_pz
                // derivative with respect to a_i
                for (int j = 0; j < 2 * degree + 1; j++) {
                    pq_d_pz(x)(i, i * (2 * degree + 2) + j) = dF(j) - dF0(j);
                }

                // derivative with respect to w
                pq_d_pz(x)(i, i * (2 * degree + 2) + (2 * degree + 1)) = kernel.dot(pdF_pw - pdF0_pw);        

                // derivative with respect to q_d0
                if (optimize_initial_velocity) {
                    if (optimize_initial_position) {
                        pq_d_pz(x)(i, Nact * (2 * degree + 2) + Nact + i) = 1;  
                    }
                    else {
                        pq_d_pz(x)(i, Nact * (2 * degree + 2) + i) = 1;
                    }
                }

                // pq_dd_pz
                // derivative with respect to a_i
                for (int j = 0; j < 2 * degree + 1; j++) {
                    pq_dd_pz(x)(i, i * (2 * degree + 2) + j) = ddF(j);
                }

                // derivative with respect to w
                pq_dd_pz(x)(i, i * (2 * degree + 2) + (2 * degree + 1)) = kernel.dot(pddF_pw);                                                                                                                                                                                            
            }
        }
    }
}

}; // namespace RAPTOR

