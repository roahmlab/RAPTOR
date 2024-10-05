#include "LMIConstraints.h"

namespace RAPTOR {

LMIConstraints::LMIConstraints(const int num_links_input,
                               const int varLength) :
    Constraints(4 * num_links_input, varLength),
    num_links(num_links_input) {
    if (10 * num_links > varLength) {
        throw std::invalid_argument("Error: regroup matrix is not provided and the variable length is too short");
    }
}

void LMIConstraints::compute(const VecX& z, 
                             bool compute_derivatives,
                             bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("Error: Hessian computation is not supported for LMI constraints");
    }        

    if (z.size() < 10 * num_links) {
        std::cerr << z.size() << std::endl;
        std::cerr << num_links << std::endl;
        throw std::invalid_argument("Error: input vector z has incorrect size");
    }

    // assume the inertial parameters are in the first 10 * num_links elements of z
    for (int i = 0; i < num_links; i++) {
        const double mass = z(10 * i);
        const Vec3& com = z.segment(10 * i + 1, 3);
        const Symmetric3 inertia(z.segment(10 * i + 4, 6));
        
        // construct LMI matrix
        Mat4 LMI = Mat4::Zero();

        LMI.topLeftCorner<3, 3>() = 
            inertia.matrix().trace() / 2.0 * Mat3::Identity() - 
            inertia.matrix();

        LMI.topRightCorner<3, 1>() = com;
        LMI.bottomLeftCorner<1, 3>() = com.transpose();

        LMI(3, 3) = mass;

        // LMI matrix semi-definiteness constraints
        // the determinant of the top-left 1x1, 2x2, 3x3, and 4x4 submatrices should be non-negative
        g(i * 4 + 0) = LMI(0, 0);
        g(i * 4 + 1) = LMI.topLeftCorner<2, 2>().determinant();
        g(i * 4 + 2) = LMI.topLeftCorner<3, 3>().determinant();
        g(i * 4 + 3) = LMI.determinant();

        if (compute_derivatives) {
            MatX pg_pz_local = MatX::Zero(4, 10);

            const double& z1 = z(10 * i + 0);
            const double& z2 = z(10 * i + 1);
            const double& z3 = z(10 * i + 2);
            const double& z4 = z(10 * i + 3);
            const double& z5 = z(10 * i + 4);
            const double& z6 = z(10 * i + 5);
            const double& z7 = z(10 * i + 6);
            const double& z8 = z(10 * i + 7);
            const double& z9 = z(10 * i + 8);
            const double& z10 = z(10 * i + 9);

            double t2 = z2*z2;
            double t3 = z3*z3;
            double t4 = z4*z4;
            double t5 = z5*z5;
            double t6 = z6*z6;
            double t7 = z7*z7;
            double t8 = z8*z8;
            double t9 = z9*z9;
            double t10 = z10*z10;
            double t11 = z2*z3*z6;
            double t12 = z2*z4*z8;
            double t13 = z3*z4*z9;
            double t14 = z5/2.0;
            double t15 = z7/2.0;
            double t16 = (z5*z7)/4.0;
            double t17 = (z5*z10)/4.0;
            double t18 = (z7*z10)/4.0;
            double t19 = -t11;
            double t20 = -t12;
            double t21 = -t13;
            double t22 = t6/2.0;
            double t23 = t8/2.0;
            double t24 = t9/2.0;
            double t25 = t16*z1;
            double t26 = t17*z1;
            double t27 = t18*z1;
            double t31 = t5/8.0;
            double t32 = t7/8.0;
            double t33 = t10/8.0;
            double t35 = t3*t14;
            double t37 = t2*t15;
            double t38 = t4*t14;
            double t40 = t4*t15;
            double t41 = (t2*z10)/2.0;
            double t42 = (t3*z10)/2.0;
            double t43 = t6*z1*(-1.0/2.0);
            double t44 = t8*z1*(-1.0/2.0);
            double t45 = t9*z1*(-1.0/2.0);
            double t28 = -t22;
            double t29 = -t23;
            double t30 = -t24;
            double t34 = t22*z1;
            double t36 = t23*z1;
            double t39 = t24*z1;
            double t46 = t31*z1;
            double t47 = t32*z1;
            double t48 = t33*z1;

            pg_pz_local(0, 4) = -1.0/2.0;
            pg_pz_local(0, 6) = 1.0/2.0;
            pg_pz_local(0, 9) = 1.0/2.0;
            pg_pz_local(1, 4) = -t14+t15;
            pg_pz_local(1, 5) = z6*-2.0;
            pg_pz_local(1, 6) = t14-t15;
            pg_pz_local(1, 9) = z10/2.0;
            pg_pz_local(2, 4) = t5*(-3.0/8.0)+t16+t17-t18+t24+t28+t29+t32+t33;
            pg_pz_local(2, 5) = -z5*z6-z6*z7+z6*z10-z8*z9*2.0;
            pg_pz_local(2, 6) = t7*(-3.0/8.0)+t16-t17+t18+t23+t28+t30+t31+t33;
            pg_pz_local(2, 7) = -z5*z8-z6*z9*2.0+z7*z8-z8*z10;
            pg_pz_local(2, 8) = z5*z9-z6*z8*2.0-z7*z9-z9*z10;
            pg_pz_local(2, 9) = t10*(-3.0/8.0)-t16+t17+t18+t22+t29+t30+t31+t32;
            pg_pz_local(3, 0) = t8*t15+t9*t14-(t6*z5)/2.0-(t6*z7)/2.0-(t8*z5)/2.0-(t9*z7)/2.0-(t8*z10)/2.0-(t9*z10)/2.0+t22*z10+t32*z5+t31*z7+t33*z5+t33*z7+t31*z10+t32*z10-(z5*z5*z5)/8.0-(z7*z7*z7)/8.0-(z10*z10*z10)/8.0-(z5*z7*z10)/4.0-z6*z8*z9*2.0;
            pg_pz_local(3, 1) = t5*z2*(-1.0/2.0)+(t7*z2)/2.0+t9*z2*2.0+(t10*z2)/2.0-z3*z5*z6-z3*z6*z7-z4*z5*z8-z2*z7*z10+z3*z6*z10-z4*z6*z9*2.0+z4*z7*z8-z3*z8*z9*2.0-z4*z8*z10;
            pg_pz_local(3, 2) = (t5*z3)/2.0-(t7*z3)/2.0+t8*z3*2.0+(t10*z3)/2.0-z2*z5*z6-z2*z6*z7+z2*z6*z10-z3*z5*z10+z4*z5*z9-z4*z6*z8*2.0-z2*z8*z9*2.0-z4*z7*z9-z4*z9*z10;
            pg_pz_local(3, 3) = (t5*z4)/2.0+t6*z4*2.0+(t7*z4)/2.0-(t10*z4)/2.0-z2*z5*z8-z4*z5*z7-z2*z6*z9*2.0+z2*z7*z8+z3*z5*z9-z3*z6*z8*2.0-z3*z7*z9-z2*z8*z10-z3*z9*z10;
            pg_pz_local(3, 4) = t13+t19+t20+t25+t26+t35+t38+t39-t42+t43+t44+t47+t48-t5*z1*(3.0/8.0)-(t2*z5)/2.0-(t4*z7)/2.0-(z1*z7*z10)/4.0;
            pg_pz_local(3, 5) = t4*z6*2.0-z2*z3*z5-z1*z5*z6-z2*z3*z7-z1*z6*z7+z2*z3*z10-z2*z4*z9*2.0-z3*z4*z8*2.0+z1*z6*z10-z1*z8*z9*2.0;
            pg_pz_local(3, 6) = t12+t19+t21+t25+t27+t36+t37+t40-t41+t43+t45+t46+t48-t7*z1*(3.0/8.0)-(t4*z5)/2.0-(t3*z7)/2.0-(z1*z5*z10)/4.0;
            pg_pz_local(3, 7) = t3*z8*2.0-z2*z4*z5+z2*z4*z7-z3*z4*z6*2.0-z1*z5*z8-z2*z3*z9*2.0-z1*z6*z9*2.0+z1*z7*z8-z2*z4*z10-z1*z8*z10;
            pg_pz_local(3, 8) = t2*z9*2.0-z2*z4*z6*2.0+z3*z4*z5-z2*z3*z8*2.0-z3*z4*z7+z1*z5*z9-z1*z6*z8*2.0-z1*z7*z9-z3*z4*z10-z1*z9*z10;
            pg_pz_local(3, 9) = t11+t20+t21+t26+t27+t34+t41+t42+t44+t45+t46+t47-(t3*z5)/2.0-(t2*z7)/2.0-t10*z1*(3.0/8.0)-(t4*z10)/2.0-(z1*z5*z7)/4.0;
        
            pg_pz.block(4 * i, 10 * i, 4, 10) = pg_pz_local;
        }
    }
}

void LMIConstraints::compute_bounds() {
    g_lb = VecX::Zero(m);
    g_ub = VecX::Constant(m, 1e19);
}

void LMIConstraints::print_violation_info() {
    for (int i = 0; i < num_links; i++) {
        for (int j = 0; j < 4; j++) {
            if (g(i * 4 + j) < 0) {
                std::cout << "LMIConstraints.cpp: determinant of submatrix " 
                          << j 
                          << " for link " 
                          << i 
                          << " is violated: " 
                          << g(i * 4 + j) 
                          << " < 0" 
                          << std::endl;
            }
        }
    }
}

}; // namespace RAPTOR