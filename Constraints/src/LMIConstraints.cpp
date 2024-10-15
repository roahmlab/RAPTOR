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
            double t43 = t6*z1*(-0.5);
            double t44 = t8*z1*(-0.5);
            double t45 = t9*z1*(-0.5);
            double t28 = -t22;
            double t29 = -t23;
            double t30 = -t24;
            double t34 = t22*z1;
            double t36 = t23*z1;
            double t39 = t24*z1;
            double t46 = t31*z1;
            double t47 = t32*z1;
            double t48 = t33*z1;

            pg_pz_local(0, 4) = -0.5;
            pg_pz_local(0, 6) = 0.5;
            pg_pz_local(0, 9) = 0.5;
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
            pg_pz_local(3, 1) = t5*z2*(-0.5)+(t7*z2)/2.0+t9*z2*2.0+(t10*z2)/2.0-z3*z5*z6-z3*z6*z7-z4*z5*z8-z2*z7*z10+z3*z6*z10-z4*z6*z9*2.0+z4*z7*z8-z3*z8*z9*2.0-z4*z8*z10;
            pg_pz_local(3, 2) = (t5*z3)/2.0-(t7*z3)/2.0+t8*z3*2.0+(t10*z3)/2.0-z2*z5*z6-z2*z6*z7+z2*z6*z10-z3*z5*z10+z4*z5*z9-z4*z6*z8*2.0-z2*z8*z9*2.0-z4*z7*z9-z4*z9*z10;
            pg_pz_local(3, 3) = (t5*z4)/2.0+t6*z4*2.0+(t7*z4)/2.0-(t10*z4)/2.0-z2*z5*z8-z4*z5*z7-z2*z6*z9*2.0+z2*z7*z8+z3*z5*z9-z3*z6*z8*2.0-z3*z7*z9-z2*z8*z10-z3*z9*z10;
            pg_pz_local(3, 4) = t13+t19+t20+t25+t26+t35+t38+t39-t42+t43+t44+t47+t48-t5*z1*(3.0/8.0)-(t2*z5)/2.0-(t4*z7)/2.0-(z1*z7*z10)/4.0;
            pg_pz_local(3, 5) = t4*z6*2.0-z2*z3*z5-z1*z5*z6-z2*z3*z7-z1*z6*z7+z2*z3*z10-z2*z4*z9*2.0-z3*z4*z8*2.0+z1*z6*z10-z1*z8*z9*2.0;
            pg_pz_local(3, 6) = t12+t19+t21+t25+t27+t36+t37+t40-t41+t43+t45+t46+t48-t7*z1*(3.0/8.0)-(t4*z5)/2.0-(t3*z7)/2.0-(z1*z5*z10)/4.0;
            pg_pz_local(3, 7) = t3*z8*2.0-z2*z4*z5+z2*z4*z7-z3*z4*z6*2.0-z1*z5*z8-z2*z3*z9*2.0-z1*z6*z9*2.0+z1*z7*z8-z2*z4*z10-z1*z8*z10;
            pg_pz_local(3, 8) = t2*z9*2.0-z2*z4*z6*2.0+z3*z4*z5-z2*z3*z8*2.0-z3*z4*z7+z1*z5*z9-z1*z6*z8*2.0-z1*z7*z9-z3*z4*z10-z1*z9*z10;
            pg_pz_local(3, 9) = t11+t20+t21+t26+t27+t34+t41+t42+t44+t45+t46+t47-(t3*z5)/2.0-(t2*z7)/2.0-t10*z1*(3.0/8.0)-(t4*z10)/2.0-(z1*z5*z7)/4.0;
        
            pg_pz.block(4 * i, 10 * i, 4, 10) = pg_pz_local;

            if (compute_hessian) {
                MatX pg1_pz_pz_local = MatX::Zero(10, 10);
                MatX pg2_pz_pz_local = MatX::Zero(10, 10);
                MatX pg3_pz_pz_local = MatX::Zero(10, 10);
                MatX pg4_pz_pz_local = MatX::Zero(10, 10);

                pg2_pz_pz_local(4, 4) = -0.5;
                pg2_pz_pz_local(4, 6) = 0.5;
                pg2_pz_pz_local(5, 5) = -2.0;
                pg2_pz_pz_local(6, 4) = 0.5;
                pg2_pz_pz_local(6, 6) = -0.5;
                pg2_pz_pz_local(9, 9) = 0.5;

                t2 = z6*2.0;
                t3 = z8*2.0;
                t4 = z9*2.0;
                t5 = -z5;
                t6 = -z6;
                t8 = -z7;
                t9 = -z8;
                t11 = -z9;
                t13 = -z10;
                t14 = z5/4.0;
                t15 = z7/4.0;
                t16 = z10/4.0;
                t7 = -t2;
                t10 = -t3;
                t12 = -t4;
                t17 = -t14;
                t18 = -t15;
                t19 = -t16;
                t20 = t15+t16+t17;
                t21 = t14+t16+t18;
                t22 = t14+t15+t19;
                pg3_pz_pz_local(4, 4) = t15+t16-z5*(3.0/4.0);
                pg3_pz_pz_local(4, 5) = t6;
                pg3_pz_pz_local(4, 6) = t22;
                pg3_pz_pz_local(4, 7) = t9;
                pg3_pz_pz_local(4, 8) = z9;
                pg3_pz_pz_local(4, 9) = t21;
                pg3_pz_pz_local(5, 4) = t6;
                pg3_pz_pz_local(5, 5) = t5+t8+z10;
                pg3_pz_pz_local(5, 6) = t6;
                pg3_pz_pz_local(5, 7) = t12;
                pg3_pz_pz_local(5, 8) = t10;
                pg3_pz_pz_local(5, 9) = z6;
                pg3_pz_pz_local(6, 4) = t22;
                pg3_pz_pz_local(6, 5) = t6;
                pg3_pz_pz_local(6, 6) = t14+t16-z7*(3.0/4.0);
                pg3_pz_pz_local(6, 7) = z8;
                pg3_pz_pz_local(6, 8) = t11;
                pg3_pz_pz_local(6, 9) = t20;
                pg3_pz_pz_local(7, 4) = t9;
                pg3_pz_pz_local(7, 5) = t12;
                pg3_pz_pz_local(7, 6) = z8;
                pg3_pz_pz_local(7, 7) = t5+t13+z7;
                pg3_pz_pz_local(7, 8) = t7;
                pg3_pz_pz_local(7, 9) = t9;
                pg3_pz_pz_local(8, 4) = z9;
                pg3_pz_pz_local(8, 5) = t10;
                pg3_pz_pz_local(8, 6) = t11;
                pg3_pz_pz_local(8, 7) = t7;
                pg3_pz_pz_local(8, 8) = t8+t13+z5;
                pg3_pz_pz_local(8, 9) = t11;
                pg3_pz_pz_local(9, 4) = t21;
                pg3_pz_pz_local(9, 5) = z6;
                pg3_pz_pz_local(9, 6) = t20;
                pg3_pz_pz_local(9, 7) = t9;
                pg3_pz_pz_local(9, 8) = t11;
                pg3_pz_pz_local(9, 9) = t14+t15-z10*(3.0/4.0);

                t2 = z2*z3;
                t3 = z1*z5;
                t4 = z2*z4;
                t5 = z1*z6;
                t6 = z2*z5;
                t7 = z3*z4;
                t8 = z1*z7;
                t9 = z2*z6;
                t10 = z3*z5;
                t11 = z1*z8;
                t12 = z2*z7;
                t13 = z3*z6;
                t14 = z4*z5;
                t15 = z1*z9;
                t16 = z2*z8;
                t17 = z3*z7;
                t18 = z1*z10;
                t19 = z4*z7;
                t20 = z5*z6;
                t21 = z2*z10;
                t22 = z3*z9;
                t23 = z4*z8;
                t24 = z3*z10;
                t25 = z4*z9;
                t26 = z5*z8;
                t27 = z6*z7;
                t28 = z4*z10;
                t29 = z5*z9;
                t30 = z7*z8;
                t31 = z6*z10;
                t32 = z7*z9;
                t33 = z8*z10;
                t34 = z9*z10;
                t35 = z2*z2;
                t36 = z3*z3;
                t37 = z4*z4;
                t38 = z5*z5;
                t39 = z6*z6;
                t40 = z7*z7;
                t41 = z8*z8;
                t42 = z9*z9;
                t43 = z10*z10;
                double t53 = z4*z6*2.0;
                double t54 = z2*z9*2.0;
                double t55 = z3*z8*2.0;
                double t57 = z4*z6*4.0;
                double t59 = z2*z9*4.0;
                double t60 = z3*z8*4.0;
                double t62 = z6*z8*2.0;
                double t63 = z6*z9*2.0;
                double t64 = z8*z9*2.0;
                double t116 = (z5*z7)/4.0;
                double t117 = (z5*z10)/4.0;
                double t118 = (z7*z10)/4.0;
                t44 = t2*2.0;
                t45 = t4*2.0;
                t46 = t5*2.0;
                t47 = t7*2.0;
                t48 = t9*2.0;
                double t49 = t11*2.0;
                double t50 = t13*2.0;
                double t51 = t15*2.0;
                double t52 = t16*2.0;
                double t56 = t22*2.0;
                double t58 = t23*2.0;
                double t61 = t25*2.0;
                double t65 = -t2;
                double t66 = -t3;
                double t68 = -t4;
                double t69 = -t5;
                double t71 = -t6;
                double t72 = -t7;
                double t74 = -t8;
                double t75 = -t9;
                double t77 = -t10;
                double t78 = -t11;
                double t80 = -t12;
                double t81 = -t13;
                double t82 = -t14;
                double t84 = -t15;
                double t85 = -t16;
                double t87 = -t17;
                double t89 = -t18;
                double t91 = -t53;
                double t92 = -t19;
                double t93 = -t20;
                double t94 = -t54;
                double t95 = -t21;
                double t96 = -t55;
                double t97 = -t22;
                double t98 = -t23;
                double t100 = -t24;
                double t102 = -t25;
                double t103 = -t26;
                double t104 = -t27;
                double t106 = -t28;
                double t107 = -t62;
                double t108 = -t63;
                double t109 = -t32;
                double t110 = -t64;
                double t111 = -t33;
                double t112 = -t34;
                double t113 = t3/4.0;
                double t114 = t8/4.0;
                double t115 = t18/4.0;
                double t119 = t2+t5;
                double t120 = t4+t11;
                double t121 = t7+t15;
                double t122 = t35/2.0;
                double t123 = t36/2.0;
                double t124 = t37/2.0;
                double t125 = t38/2.0;
                double t126 = t39/2.0;
                double t127 = t40/2.0;
                double t128 = t41/2.0;
                double t129 = t42/2.0;
                double t130 = t43/2.0;
                double t134 = -t116;
                double t135 = -t117;
                double t136 = -t118;
                double t143 = t38/8.0;
                double t144 = t38*(3.0/8.0);
                double t145 = t40/8.0;
                double t146 = t40*(3.0/8.0);
                double t147 = t43/8.0;
                double t148 = t43*(3.0/8.0);
                double t67 = -t44;
                double t70 = -t45;
                double t73 = -t46;
                double t76 = -t47;
                double t79 = -t48;
                double t83 = -t49;
                double t86 = -t50;
                double t88 = -t51;
                double t90 = -t52;
                double t99 = -t56;
                double t101 = -t58;
                double t105 = -t61;
                double t131 = -t113;
                double t132 = -t114;
                double t133 = -t115;
                double t137 = -t122;
                double t138 = -t123;
                double t139 = -t124;
                double t140 = -t126;
                double t141 = -t128;
                double t142 = -t129;
                double t149 = -t144;
                double t150 = -t146;
                double t151 = -t148;
                double t152 = t65+t69;
                double t154 = t68+t78;
                double t156 = t72+t84;
                double t158 = t60+t91+t94;
                double t159 = t59+t91+t96;
                double t160 = t57+t94+t96;
                double t161 = t71+t81+t98;
                double t162 = t75+t87+t102;
                double t163 = t85+t97+t106;
                double t164 = t9+t24+t77+t102;
                double t165 = t12+t23+t81+t95;
                double t166 = t14+t22+t85+t92;
                double t167 = t10+t25+t75+t100;
                double t168 = t13+t21+t80+t98;
                double t169 = t16+t19+t82+t97;
                double t176 = t31+t93+t104+t110;
                double t177 = t30+t103+t108+t111;
                double t178 = t29+t107+t109+t112;
                double t153 = t67+t73;
                double t155 = t70+t83;
                double t157 = t76+t88;
                double t170 = t12+t71+t86+t95;
                double t171 = t10+t79+t87+t100;
                double t172 = t21+t71+t80+t101;
                double t173 = t14+t90+t92+t106;
                double t174 = t24+t77+t87+t105;
                double t175 = t19+t82+t99+t106;
                double t179 = t114+t115+t131+t137;
                double t180 = t113+t115+t132+t138;
                double t181 = t113+t114+t133+t139;
                double t182 = t117+t118+t126+t134+t141+t142+t143+t145+t151;
                double t183 = t116+t118+t128+t135+t140+t142+t143+t147+t150;
                double t184 = t116+t117+t129+t136+t140+t141+t145+t147+t149;
                pg4_pz_pz_local(0, 4) = t184;
                pg4_pz_pz_local(0, 5) = t176;
                pg4_pz_pz_local(0, 6) = t183;
                pg4_pz_pz_local(0, 7) = t177;
                pg4_pz_pz_local(0, 8) = t178;
                pg4_pz_pz_local(0, 9) = t182;
                pg4_pz_pz_local(1, 1) = t42*2.0-t125+t127+t130-z7*z10;
                pg4_pz_pz_local(1, 2) = t176;
                pg4_pz_pz_local(1, 3) = t177;
                pg4_pz_pz_local(1, 4) = t161;
                pg4_pz_pz_local(1, 5) = t174;
                pg4_pz_pz_local(1, 6) = t165;
                pg4_pz_pz_local(1, 7) = t175;
                pg4_pz_pz_local(1, 8) = t159;
                pg4_pz_pz_local(1, 9) = t168;
                pg4_pz_pz_local(2, 1) = t176;
                pg4_pz_pz_local(2, 2) = t41*2.0+t125-t127+t130-z5*z10;
                pg4_pz_pz_local(2, 3) = t178;
                pg4_pz_pz_local(2, 4) = t167;
                pg4_pz_pz_local(2, 5) = t172;
                pg4_pz_pz_local(2, 6) = t162;
                pg4_pz_pz_local(2, 7) = t158;
                pg4_pz_pz_local(2, 8) = t173;
                pg4_pz_pz_local(2, 9) = t164;
                pg4_pz_pz_local(3, 1) = t177;
                pg4_pz_pz_local(3, 2) = t178;
                pg4_pz_pz_local(3, 3) = t39*2.0+t125+t127-t130-z5*z7;
                pg4_pz_pz_local(3, 4) = t166;
                pg4_pz_pz_local(3, 5) = t160;
                pg4_pz_pz_local(3, 6) = t169;
                pg4_pz_pz_local(3, 7) = t170;
                pg4_pz_pz_local(3, 8) = t171;
                pg4_pz_pz_local(3, 9) = t163;
                pg4_pz_pz_local(4, 0) = t184;
                pg4_pz_pz_local(4, 1) = t161;
                pg4_pz_pz_local(4, 2) = t167;
                pg4_pz_pz_local(4, 3) = t166;
                pg4_pz_pz_local(4, 4) = t3*(-3.0/4.0)+t114+t115+t123+t124+t137;
                pg4_pz_pz_local(4, 5) = t152;
                pg4_pz_pz_local(4, 6) = t181;
                pg4_pz_pz_local(4, 7) = t154;
                pg4_pz_pz_local(4, 8) = t121;
                pg4_pz_pz_local(4, 9) = t180;
                pg4_pz_pz_local(5, 0) = t176;
                pg4_pz_pz_local(5, 1) = t174;
                pg4_pz_pz_local(5, 2) = t172;
                pg4_pz_pz_local(5, 3) = t160;
                pg4_pz_pz_local(5, 4) = t152;
                pg4_pz_pz_local(5, 5) = t18+t37*2.0+t66+t74;
                pg4_pz_pz_local(5, 6) = t152;
                pg4_pz_pz_local(5, 7) = t157;
                pg4_pz_pz_local(5, 8) = t155;
                pg4_pz_pz_local(5, 9) = t119;
                pg4_pz_pz_local(6, 0) = t183;
                pg4_pz_pz_local(6, 1) = t165;
                pg4_pz_pz_local(6, 2) = t162;
                pg4_pz_pz_local(6, 3) = t169;
                pg4_pz_pz_local(6, 4) = t181;
                pg4_pz_pz_local(6, 5) = t152;
                pg4_pz_pz_local(6, 6) = t8*(-3.0/4.0)+t113+t115+t122+t124+t138;
                pg4_pz_pz_local(6, 7) = t120;
                pg4_pz_pz_local(6, 8) = t156;
                pg4_pz_pz_local(6, 9) = t179;
                pg4_pz_pz_local(7, 0) = t177;
                pg4_pz_pz_local(7, 1) = t175;
                pg4_pz_pz_local(7, 2) = t158;
                pg4_pz_pz_local(7, 3) = t170;
                pg4_pz_pz_local(7, 4) = t154;
                pg4_pz_pz_local(7, 5) = t157;
                pg4_pz_pz_local(7, 6) = t120;
                pg4_pz_pz_local(7, 7) = t8+t36*2.0+t66+t89;
                pg4_pz_pz_local(7, 8) = t153;
                pg4_pz_pz_local(7, 9) = t154;
                pg4_pz_pz_local(8, 0) = t178;
                pg4_pz_pz_local(8, 1) = t159;
                pg4_pz_pz_local(8, 2) = t173;
                pg4_pz_pz_local(8, 3) = t171;
                pg4_pz_pz_local(8, 4) = t121;
                pg4_pz_pz_local(8, 5) = t155;
                pg4_pz_pz_local(8, 6) = t156;
                pg4_pz_pz_local(8, 7) = t153;
                pg4_pz_pz_local(8, 8) = t3+t35*2.0+t74+t89;
                pg4_pz_pz_local(8, 9) = t156;
                pg4_pz_pz_local(9, 0) = t182;
                pg4_pz_pz_local(9, 1) = t168;
                pg4_pz_pz_local(9, 2) = t164;
                pg4_pz_pz_local(9, 3) = t163;
                pg4_pz_pz_local(9, 4) = t180;
                pg4_pz_pz_local(9, 5) = t119;
                pg4_pz_pz_local(9, 6) = t179;
                pg4_pz_pz_local(9, 7) = t154;
                pg4_pz_pz_local(9, 8) = t156;
                pg4_pz_pz_local(9, 9) = t18*(-3.0/4.0)+t113+t114+t122+t123+t139;

                pg_pz_pz(4 * i + 0).block(10 * i, 10 * i, 10, 10) = pg1_pz_pz_local;
                pg_pz_pz(4 * i + 1).block(10 * i, 10 * i, 10, 10) = pg2_pz_pz_local;
                pg_pz_pz(4 * i + 2).block(10 * i, 10 * i, 10, 10) = pg3_pz_pz_local;
                pg_pz_pz(4 * i + 3).block(10 * i, 10 * i, 10, 10) = pg4_pz_pz_local;
            }
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