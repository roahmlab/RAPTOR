#include "LMIConstraints.h"

namespace RAPTOR {

LMIConstraints::LMIConstraints(const int num_links_input,
                               const int varLength) :
    Constraints(10 * num_links_input, varLength),
    num_links(num_links_input) {
    if (10 * num_links > varLength) {
        throw std::invalid_argument("Error: regroup matrix is not provided and the variable length is too short");
    }

    // initialize d LMI d z
    dLMIdz(0) = MatX::Zero(4, 4);
    dLMIdz(0)(3, 3) = 1;

    dLMIdz(1) = MatX::Zero(4, 4);
    dLMIdz(1)(0, 3) = 1;
    dLMIdz(1)(3, 0) = 1;

    dLMIdz(2) = MatX::Zero(4, 4);
    dLMIdz(2)(1, 3) = 1;
    dLMIdz(2)(3, 1) = 1;

    dLMIdz(3) = MatX::Zero(4, 4);
    dLMIdz(3)(2, 3) = 1;
    dLMIdz(3)(3, 2) = 1;

    dLMIdz(4) = MatX::Zero(4, 4);
    dLMIdz(4)(0, 0) = -0.5;
    dLMIdz(4)(1, 1) = 0.5;
    dLMIdz(4)(2, 2) = 0.5;

    dLMIdz(5) = MatX::Zero(4, 4);
    dLMIdz(5)(0, 1) = -1;
    dLMIdz(5)(1, 0) = -1;

    dLMIdz(6) = MatX::Zero(4, 4);
    dLMIdz(6)(0, 0) = 0.5;
    dLMIdz(6)(1, 1) = -0.5;
    dLMIdz(6)(2, 2) = -0.5;

    dLMIdz(7) = MatX::Zero(4, 4);
    dLMIdz(7)(0, 2) = -1;
    dLMIdz(7)(2, 0) = -1;

    dLMIdz(8) = MatX::Zero(4, 4);   
    dLMIdz(8)(1, 2) = -1;
    dLMIdz(8)(2, 1) = -1;

    dLMIdz(9) = MatX::Zero(4, 4);
    dLMIdz(9)(0, 0) = 0.5;
    dLMIdz(9)(1, 1) = 0.5;
    dLMIdz(9)(2, 2) = -0.5;
}

void LMIConstraints::compute(const VecX& z, 
                             bool compute_derivatives,
                             bool compute_hessian) {  
    if (compute_hessian) {
        throw std::invalid_argument("Error: LMIConstraints does not support hessian computation");
    }

    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
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
            0.5 * inertia.matrix().trace() * Mat3::Identity() - 
            inertia.matrix();

        LMI.topRightCorner<3, 1>() = com;
        LMI.bottomLeftCorner<1, 3>() = com.transpose();

        LMI(3, 3) = mass;

        // LMI matrix semi-definiteness constraints
        // the principal minors of the LMI matrix should be non-negative
        g(i * 10 + 0) = LMI(0, 0);
        g(i * 10 + 1) = LMI(1, 1);
        g(i * 10 + 2) = LMI(2, 2);
        g(i * 10 + 3) = LMI(3, 3);
        g(i * 10 + 4) = LMI.topLeftCorner<2, 2>().determinant();
        g(i * 10 + 5) = LMI.block(1, 1, 2, 2).determinant();
        g(i * 10 + 6) = LMI.bottomRightCorner<2, 2>().determinant();
        g(i * 10 + 7) = LMI.topLeftCorner<3, 3>().determinant();
        g(i * 10 + 8) = LMI.bottomRightCorner<3, 3>().determinant();
        g(i * 10 + 9) = LMI.determinant();

        if (compute_derivatives) {MatX pg_pz_local = MatX::Zero(10, 10);

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
            double t14 = z1/2.0;
            double t15 = z5/2.0;
            double t16 = z7/2.0;
            double t17 = z10/2.0;
            double t20 = (z5*z7)/4.0;
            double t21 = (z5*z10)/4.0;
            double t22 = (z7*z10)/4.0;
            double t18 = t14*z7;
            double t19 = t14*z10;
            double t23 = -t11;
            double t24 = -t12;
            double t25 = -t13;
            double t26 = -t16;
            double t27 = -t17;
            double t28 = t3/2.0;
            double t29 = t4/2.0;
            double t30 = t6/2.0;
            double t31 = t8/2.0;
            double t32 = t9/2.0;
            double t33 = t20*z1;
            double t34 = t21*z1;
            double t35 = t22*z1;
            double t41 = t5/8.0;
            double t42 = t7/8.0;
            double t43 = t10/8.0;
            double t44 = t6*t14;
            double t45 = t3*t15;
            double t46 = t8*t14;
            double t47 = t2*t16;
            double t48 = t4*t15;
            double t49 = t9*t14;
            double t50 = t4*t16;
            double t51 = t2*t17;
            double t52 = t3*t17;
            double t53 = t6*z1*(-1.0/2.0);
            double t54 = t8*z1*(-1.0/2.0);
            double t55 = t9*z1*(-1.0/2.0);
            double t36 = -t28;
            double t37 = -t29;
            double t38 = -t30;
            double t39 = -t31;
            double t40 = -t32;
            double t56 = t41*z1;
            double t57 = t42*z1;
            double t58 = t43*z1;
            pg_pz_local(0, 4) = -1.0/2.0;
            pg_pz_local(0, 6) = 1.0/2.0;
            pg_pz_local(0, 9) = 1.0/2.0;
            pg_pz_local(1, 4) = 1.0/2.0;
            pg_pz_local(1, 6) = -1.0/2.0;
            pg_pz_local(1, 9) = 1.0/2.0;
            pg_pz_local(2, 4) = 1.0/2.0;
            pg_pz_local(2, 6) = 1.0/2.0;
            pg_pz_local(2, 9) = -1.0/2.0;
            pg_pz_local(3, 0) = 1.0;
            pg_pz_local(4, 4) = -t15+t16;
            pg_pz_local(4, 5) = z6*-2.0;
            pg_pz_local(4, 6) = t15+t26;
            pg_pz_local(4, 9) = t17;
            pg_pz_local(5, 4) = t15;
            pg_pz_local(5, 6) = t17+t26;
            pg_pz_local(5, 8) = z9*-2.0;
            pg_pz_local(5, 9) = t16+t27;
            pg_pz_local(6, 0) = t15+t16+t27;
            pg_pz_local(6, 3) = z4*-2.0;
            pg_pz_local(6, 4) = t14;
            pg_pz_local(6, 6) = t14;
            pg_pz_local(6, 9) = -t14;
            pg_pz_local(7, 4) = t5*(-3.0/8.0)+t20+t21-t22+t32+t38+t39+t42+t43;
            pg_pz_local(7, 5) = -z5*z6-z6*z7+z6*z10-z8*z9*2.0;
            pg_pz_local(7, 6) = t7*(-3.0/8.0)+t20-t21+t22+t31+t38+t40+t41+t43;
            pg_pz_local(7, 7) = -z5*z8-z6*z9*2.0+z7*z8-z8*z10;
            pg_pz_local(7, 8) = z5*z9-z6*z8*2.0-z7*z9-z9*z10;
            pg_pz_local(7, 9) = t10*(-3.0/8.0)-t20+t21+t22+t30+t39+t40+t41+t42;
            pg_pz_local(8, 0) = t5/4.0-t7/4.0-t9-t10/4.0+t16*z10;
            pg_pz_local(8, 2) = -z3*z5-z3*z7+z3*z10-z4*z9*2.0;
            pg_pz_local(8, 3) = -z4*z5+z4*z7-z3*z9*2.0-z4*z10;
            pg_pz_local(8, 4) = t36+t37+t14*z5;
            pg_pz_local(8, 6) = t19+t29+t36-(z1*z7)/2.0;
            pg_pz_local(8, 8) = z3*z4*-2.0-z1*z9*2.0;
            pg_pz_local(8, 9) = t18+t28+t37-(z1*z10)/2.0;
            pg_pz_local(9, 0) = t6*t17+t8*t16+t9*t15-(t6*z5)/2.0-(t6*z7)/2.0-(t8*z5)/2.0-(t9*z7)/2.0-(t8*z10)/2.0-(t9*z10)/2.0+t42*z5+t41*z7+t43*z5+t43*z7+t41*z10+t42*z10-(z5*z5*z5)/8.0-(z7*z7*z7)/8.0-(z10*z10*z10)/8.0-(z5*z7*z10)/4.0-z6*z8*z9*2.0;
            pg_pz_local(9, 1) = t5*z2*(-1.0/2.0)+(t7*z2)/2.0+t9*z2*2.0+(t10*z2)/2.0-z3*z5*z6-z3*z6*z7-z4*z5*z8-z2*z7*z10+z3*z6*z10-z4*z6*z9*2.0+z4*z7*z8-z3*z8*z9*2.0-z4*z8*z10;
            pg_pz_local(9, 2) = (t5*z3)/2.0-(t7*z3)/2.0+t8*z3*2.0+(t10*z3)/2.0-z2*z5*z6-z2*z6*z7+z2*z6*z10-z3*z5*z10+z4*z5*z9-z4*z6*z8*2.0-z2*z8*z9*2.0-z4*z7*z9-z4*z9*z10;
            pg_pz_local(9, 3) = (t5*z4)/2.0+t6*z4*2.0+(t7*z4)/2.0-(t10*z4)/2.0-z2*z5*z8-z4*z5*z7-z2*z6*z9*2.0+z2*z7*z8+z3*z5*z9-z3*z6*z8*2.0-z3*z7*z9-z2*z8*z10-z3*z9*z10;
            pg_pz_local(9, 4) = t13+t23+t24+t33+t34+t45+t48+t49+t53+t54+t57+t58-t5*z1*(3.0/8.0)-(t2*z5)/2.0-(t4*z7)/2.0-(t3*z10)/2.0-(z1*z7*z10)/4.0;
            pg_pz_local(9, 5) = t4*z6*2.0-z2*z3*z5-z1*z5*z6-z2*z3*z7-z1*z6*z7+z2*z3*z10-z2*z4*z9*2.0-z3*z4*z8*2.0+z1*z6*z10-z1*z8*z9*2.0;
            pg_pz_local(9, 6) = t12+t23+t25+t33+t35+t46+t47+t50+t53+t55+t56+t58-t7*z1*(3.0/8.0)-(t4*z5)/2.0-(t3*z7)/2.0-(t2*z10)/2.0-(z1*z5*z10)/4.0;
            pg_pz_local(9, 7) = t3*z8*2.0-z2*z4*z5+z2*z4*z7-z3*z4*z6*2.0-z1*z5*z8-z2*z3*z9*2.0-z1*z6*z9*2.0+z1*z7*z8-z2*z4*z10-z1*z8*z10;
            pg_pz_local(9, 8) = t2*z9*2.0-z2*z4*z6*2.0+z3*z4*z5-z2*z3*z8*2.0-z3*z4*z7+z1*z5*z9-z1*z6*z8*2.0-z1*z7*z9-z3*z4*z10-z1*z9*z10;
            pg_pz_local(9, 9) = t11+t24+t25+t34+t35+t44+t51+t52+t54+t55+t56+t57-(t3*z5)/2.0-(t2*z7)/2.0-t10*z1*(3.0/8.0)-(t4*z10)/2.0-(z1*z5*z7)/4.0;

            pg_pz.block(10 * i, 10 * i, 10, 10) = pg_pz_local;
        }
    }
}

void LMIConstraints::compute_bounds() {
    g_lb = VecX::Constant(m, 1e-3);
    g_ub = VecX::Constant(m, 1e19);
}

void LMIConstraints::print_violation_info() {
    for (int i = 0; i < num_links; i++) {
        for (int j = 0; j < 10; j++) {
            if (g(i * 10 + j) < 0) {
                std::cout << "LMIConstraints.cpp: determinant of submatrix " 
                          << j 
                          << " for link " 
                          << i 
                          << " is violated: " 
                          << g(i * 10 + j) 
                          << " < 0" 
                          << std::endl;
            }
        }
    }
}

}; // namespace RAPTOR