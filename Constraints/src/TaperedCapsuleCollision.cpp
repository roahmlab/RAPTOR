#include "TaperedCapsuleCollision.h"

namespace RAPTOR {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;

double solve_quadratic(double a, double b, double c, int sign){
    return (-b+sign*sqrt(pow(b,2)-4*a*c))/(2*a);
}

MatX batchDot(const Vec3 vector, const MatX matrix){
    int dims = matrix.cols();
    VecX result(dims);
    
    for(uint i = 0; i<dims; i++){
        // std::cout << "Vector\n" << vector << "Matrix\n" << matrix.col(i) << "\n";
        result[i] = vector.dot(matrix.col(i));
    }
    return result;
}

double distanceInternal(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2,
                        MatX& pdist_pz){
    

    double r1 = tc1_radius_2-tc1_radius_1;
    double r2 = tc2_radius_2-tc2_radius_1;

    Vec3 d1 = tc1_point_2-tc1_point_1;
    Vec3 d2 = tc2_point_2-tc2_point_1;
    Vec3 d12 = tc2_point_1-tc1_point_1;

    double a = d1.dot(d1);
    double b = d1.dot(d2);
    double c = d1.dot(d12);
    double e = d2.dot(d2);
    double f = d2.dot(d12);
    double g = d12.dot(d12);

    VecX distances(14);
    VecX u_test(14);
    VecX t_test(14);
    
    // end point to end point checks
    u_test[10] = 0.0;
    u_test[11] = 0.0;
    u_test[12] = 1.0;
    u_test[13] = 1.0;
    t_test[10] = 0.0;
    t_test[11] = 1.0;
    t_test[12] = 0.0;
    t_test[13] = 1.0;
    
    // end point to edge checks
        // case 1-2
    u_test[0] = solve_quadratic(e*(e-r2*r2), 2*f*(e-r2*r2), f*f-r2*r2*g, 1);
    u_test[1] = solve_quadratic(e*(e-r2*r2), 2*f*(e-r2*r2), f*f-r2*r2*g, -1);
    u_test[2] = solve_quadratic(e*(e-r2*r2), 2*(f-b)*(e-r2*r2), (f-b)*(f-b)-r2*r2*(g-2*c+a), 1);
    u_test[3] = solve_quadratic(e*(e-r2*r2), 2*(f-b)*(e-r2*r2), (f-b)*(f-b)-r2*r2*(g-2*c+a), -1);
    t_test[0] = 0;
    t_test[1] = 0;
    t_test[2] = 1;
    t_test[3] = 1;
        // case 3-4
    u_test[4] = 0;
    u_test[5] = 0;
    u_test[6] = 1;
    u_test[7] = 1;
    t_test[4]  = solve_quadratic(a*(a-pow(r1,2)), -2*c*(a-pow(r1,2)), pow(c,2)-pow(r1,2)*g, 1);
    t_test[5]  = solve_quadratic(a*(a-pow(r1,2)), -2*c*(a-pow(r1,2)), pow(c,2)-pow(r1,2)*g, -1);
    t_test[6] = solve_quadratic(a*(a-pow(r1,2)), -2*(b+c)*(a-pow(r1,2)), pow((b+c),2)-pow(r1,2)*(e+2*f+g), 1);
    t_test[7] = solve_quadratic(a*(a-pow(r1,2)), -2*(b+c)*(a-pow(r1,2)), pow((b+c),2)-pow(r1,2)*(e+2*f+g), -1);
    // edge to edge check
        // case 5
    double phi = (r2*b+r1*f)/(r2*a+r1*b);
    double gamma = (r2*c+r1*f)/(r2*a+r1*b);
    double alpha = a*phi-b;
    double beta = a*gamma-c;
    
    double a_5 = alpha*alpha - r1*r1*phi*phi*a - r1*r1*e + 2*r1*r1*b*phi;
    double b_5 = 2*alpha*beta - 2*r1*r1*phi*gamma*a + 2*r1*r1*b*gamma + 2*r1*c*phi - 2*r1*r1*f;
    double c_5 = beta*beta - r1*r1*a*gamma*gamma + 2*r1*r1*c*gamma - r1*r1*g;
    u_test[8] = solve_quadratic(a_5, b_5, c_5, 1);
    u_test[9] = solve_quadratic(a_5, b_5, c_5, -1);
    t_test[8] = phi*u_test[8]+gamma;
    t_test[9] = phi*u_test[9]+gamma;

    t_test[8] = (c*r2 + f*r1)/(a*r2 + b*r1) + ((b*r2 + e*r1)*(a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*b*b*c + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));
    t_test[9] = (c*r2 + f*r1)/(a*r2 + b*r1) - ((b*r2 + e*r1)*(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));

    u_test[8] = (a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*b*b*c + b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);
    u_test[9] = -(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);

    for(int i = 0; i<14; i++){
        if(u_test[i] <= 1 && u_test[i] >= 0 && t_test[i] <= 1 && t_test[i] >= 0){
            distances[i] = (d1*t_test[i]-d2*u_test[i]-d12).norm() - r1*t_test[i] - r2*u_test[i];
        }
        else{
            distances[i] = INFINITY;
        }
    }
    Eigen::VectorXd::Index ind = -1;
    double finalDistance = distances.minCoeff(&ind) - tc1_radius_1 - tc2_radius_1;
    double uStar = u_test[ind];
    double tStar = t_test[ind];

    if(ptc1_point_1_pz.size()>0){
        // Solve for gradients of T and U
        MatX pd1_pz = ptc1_point_2_pz - ptc1_point_1_pz;
        MatX pd2_pz = ptc2_point_2_pz - ptc2_point_1_pz;
        MatX pd12_pz = ptc2_point_1_pz - ptc1_point_1_pz;

        std::cout << "Derivative" << "\n";
        MatX pa_pz = batchDot(d1,pd1_pz)+batchDot(d1, pd1_pz);
        MatX pb_pz = batchDot(d1,pd2_pz)+batchDot(d2, pd1_pz);
        MatX pc_pz = batchDot(d12,pd1_pz)+batchDot(d1, pd12_pz);
        MatX pe_pz = batchDot(d2,pd2_pz)+batchDot(d2, pd2_pz);
        MatX pf_pz = batchDot(d2,pd12_pz)+batchDot(d12, pd2_pz);
        MatX pg_pz = batchDot(d12,pd12_pz)+batchDot(d12, pd12_pz);

        std::cout << "pa_pz: " << pa_pz << "\n";
        std::cout << "pb_pz: " << pb_pz << "\n";
        std::cout << "pc_pz: " << pc_pz << "\n";
        std::cout << "pe_pz: " << pe_pz << "\n";
        std::cout << "pf_pz: " << pf_pz << "\n";
        std::cout << "pg_pz: " << pg_pz << "\n";

        MatX lgrey = pa_pz * tStar - pb_pz * uStar - pc_pz;
        double purple = (d1 * tStar - d2 * uStar - d12).norm();
        double orange = a * tStar - b * uStar - c;
        double green = -b * tStar - e * uStar + f;
        MatX red = -pb_pz * tStar + pe_pz * uStar + pf_pz;
        MatX grey = (batchDot((d1 * tStar*tStar - d2 * tStar * uStar - d12 * tStar),pd1_pz)-
                    batchDot((d1 * tStar * uStar - d2 * uStar*uStar - d12 * uStar),pd2_pz)-
                    batchDot((d1 - d2 - d12), pd12_pz))/ (purple*purple);

        double L11 = a / purple - orange*orange/pow(purple,3);
        double L22 = e / purple - green*green/pow(purple,3);
        double L12 = -b / purple - orange * green /pow(purple,3);

        std::cout << "L11: " << L11 << "\n";
        std::cout << "L22: " << L22 << "\n";
        std::cout << "L12: " << L12 << "\n";

        MatX dLambda1 = (tStar == 1) * -(lgrey * purple - orange * grey);
        MatX dLambda2 = (uStar == 1) * -(red * purple - green * grey);
        MatX dLambda3 = (tStar == 0) * (lgrey * purple - orange * grey);
        MatX dLambda4 = (uStar == 0) * (red * purple - green * grey);

        std::cout << "dLambda1: " << dLambda1 << "\n";

        MatX H1 = (
            (pa_pz * tStar - pb_pz * uStar - pc_pz) / purple + dLambda1 - dLambda3
            - orange * batchDot(-d1 * tStar + d2 * uStar + d12,-pd1_pz * tStar + pd2_pz * uStar + pd12_pz) / pow(purple,3)
        ).transpose();
        MatX H2 = (
            (-pb_pz * tStar + pe_pz * uStar + pf_pz) / purple  + dLambda2 - dLambda4
             - green * batchDot(-d1 * tStar + d2 * uStar + d12,-pd1_pz * tStar + pd2_pz * uStar + pd12_pz) / pow(purple,3)
        ).transpose();

        std::cout << "H1: " << H1 << "\n";
        std::cout << "H2: " << H2 << "\n";

        MatX pt_pz = MatX::Zero(1, ptc1_point_1_pz.cols());
        MatX pu_pz = MatX::Zero(1, ptc1_point_1_pz.cols());
        switch(ind){
            case 0:
            case 1:
            case 2:
            case 3:
                pt_pz = -H1 / L11;
                // pu_pz = 0;
                break;
            case 4:
            case 5:
            case 6:
            case 7:
                // pt_pz = 0;
                pu_pz = -H2 / L22;
                break;
            case 8:
            case 9:
                pt_pz = (L12 * H2 - L22 * H1) / (L11 * L22 - L12*L12);
                pu_pz = (L12 * H2 - L22 * H1) / (L11 * L22 - L12*L12);
                break;
            case 10:
            case 11:
            case 12:
            case 13:
                // pt_pz = 0;
                // pu_pz = 0;
                break;
        }

        

        std::cout << "pt_pz: " << pt_pz << "\n";
        std::cout << "pu_pz: " << pt_pz << "\n";
        
        MatX part1 = batchDot((d1 * tStar - d2 * uStar - d12), pd1_pz*tStar + d1*pt_pz - pd2_pz*uStar - d2*pu_pz - pd12_pz);
        std::cout << "Part1: " << part1 << "\n";
        pdist_pz = part1.transpose() / purple - r1 * pt_pz - r2 * pu_pz;

    }

    return finalDistance;
}

double TaperedCapsuleCollision::computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2,
                        MatX& pdist_pz){
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                ptc1_point_1_pz, ptc1_point_2_pz, ptc2_point_1_pz, ptc2_point_2_pz,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, pdist_pz);
    return distance;
}

double TaperedCapsuleCollision::computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const double tc1_radius_1, const double tc1_radius_2, 
                        const double tc2_radius_1, const double tc2_radius_2){
    MatX nullDerivative;
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                nullDerivative, nullDerivative, nullDerivative, nullDerivative,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, nullDerivative);
    return distance;
}


} // namespace RAPTOR