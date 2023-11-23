#include "FourierCurves.h"

namespace IDTO {

FourierCurves::FourierCurves(const VecX& tspan_input, int Nact_input, int degree_input) : 
    Trajectories(tspan_input, Nact_input),
    degree(degree_input) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd dF = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd ddF = Eigen::VectorXd::Zero(2 * degree + 1);

    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd dF0 = Eigen::VectorXd::Zero(2 * degree + 1);
}

FourierCurves::FourierCurves(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input) :
    Trajectories(T_input, N_input, Nact_input, time_discretization),
    degree(degree_input) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd dF = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd ddF = Eigen::VectorXd::Zero(2 * degree + 1);

    Eigen::VectorXd F0 = Eigen::VectorXd::Zero(2 * degree + 1);
    Eigen::VectorXd dF0 = Eigen::VectorXd::Zero(2 * degree + 1);
}

void FourierCurves::compute(const VecX& z, bool compute_derivatives) {
    Eigen::MatrixXd temp = z.head((2 * degree + 2) * Nact);
    MatX coefficients = temp.reshaped(2 * degree + 2, Nact);
    VecX q_act0       = z.block((2 * degree + 2) * Nact, 0, Nact, 1);
    VecX q_act_d0     = z.block((2 * degree + 2) * Nact + Nact, 0, Nact, 1);

    for (int x = 0; x < N; x++) {
        double t = tspan(x);

        for (int i = 0; i < Nact; i++) {
            Eigen::VectorXd kernel = coefficients.block(0, i, 2 * degree + 1, 1);
            double w = coefficients(2 * degree + 1, i);

            ddF(0) = 1;
            dF(0)  = t;
            F(0)   = t * t * 0.5;

            for (int j = 0; j < degree; j++) {
                double sinjwt = sin((j + 1) * w * t);
                double cosjwt = cos((j + 1) * w * t);

                ddF(2 * j + 1) = cosjwt;
                ddF(2 * j + 2) = sinjwt;

                double jw = (j + 1) * w;
                dF(2 * j + 1) = sinjwt / jw;
                dF(2 * j + 2) = -cosjwt / jw;
                dF0(2 * j + 2) = -1 / jw;

                double jw2 = jw * jw;
                F(2 * j + 1) = -cosjwt / jw2;
                F(2 * j + 2) = -sinjwt / jw2;
                F0(2 * j + 1) = -1 / jw2;
            }

            q_dd(x)(i) = ddF.dot(kernel);

            double q_d_raw = dF.dot(kernel);
            double q_d_raw0 = dF0.dot(kernel);
            q_d(x)(i) = q_d_raw + (q_act_d0(i) - q_d_raw0);

            double q_raw = F.dot(kernel) + (q_act_d0(i) - q_d_raw0) * t;
            double q_raw0 = F0.dot(kernel);
            q(x)(i) = q_raw + (q_act0(i) - q_raw0);
        }
    }
}

}; // namespace IDTO

