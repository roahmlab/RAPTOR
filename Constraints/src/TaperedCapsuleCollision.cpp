#include "TaperedCapsuleCollision.h"
#include <chrono>

namespace RAPTOR {

inline double solve_quadratic(double a, double b, double c, int sign){
    return (-b+sign*sqrt(pow(b,2)-4*a*c))/(2*a);
}

inline Eigen::Vector<double,NUM_FACTORS> batchDot(const Eigen::Vector3d vector, const Eigen::Matrix<double,3,NUM_FACTORS> matrix){
    // equivalent to this line, but for loop is faster
    // (matrix.transpose()*vector)
    int dims = matrix.cols();
    Eigen::VectorXd result(dims);
    
    for(uint i = 0; i<dims; i++){
        result[i] = vector.dot(matrix.col(i));
    }
    return result;
}

inline double TaperedCapsuleCollision::distanceInternal(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double& tc1_radius_1, const double& tc1_radius_2, 
                        const double& tc2_radius_1, const double& tc2_radius_2,
                        VecX& pdist_pz){
    // auto start = std::chrono::high_resolution_clock::now();
    // std::cout << "Entering func" << std::endl;

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

    Eigen::Vector<double, 14> distances(14);
    Eigen::Vector<double, 14> u_test(14);
    Eigen::Vector<double, 14> t_test(14);
    
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
    // double phi = (r2*b+r1*f)/(r2*a+r1*b);
    // double gamma = (r2*c+r1*f)/(r2*a+r1*b);
    // double alpha = a*phi-b;
    // double beta = a*gamma-c;
    
    // double a_5 = alpha*alpha - r1*r1*phi*phi*a - r1*r1*e + 2*r1*r1*b*phi;
    // double b_5 = 2*alpha*beta - 2*r1*r1*phi*gamma*a + 2*r1*r1*b*gamma + 2*r1*c*phi - 2*r1*r1*f;
    // double c_5 = beta*beta - r1*r1*a*gamma*gamma + 2*r1*r1*c*gamma - r1*r1*g;
    // u_test[8] = solve_quadratic(a_5, b_5, c_5, 1);
    // u_test[9] = solve_quadratic(a_5, b_5, c_5, -1);
    // t_test[8] = phi*u_test[8]+gamma;
    // t_test[9] = phi*u_test[9]+gamma;

    

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

    // std::cout << "Finding dist" << std::endl;

    Eigen::VectorXd::Index ind = -1;
    double finalDistance = distances.minCoeff(&ind) - tc1_radius_1 - tc2_radius_1;
    double uStar = u_test[ind];
    double tStar = t_test[ind];

    // auto time1 = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time1 - start);
    // std::cout << "Time taken to find dist: " << duration.count()/1e3 << " microseconds" << std::endl;

    if(ptc1_point_1_pz.size()>0){
        int colmns = ptc1_point_1_pz.cols();
        // Solve for gradients of T and U
        MatX pd1_pz = ptc1_point_2_pz - ptc1_point_1_pz;
        MatX pd2_pz = ptc2_point_2_pz - ptc2_point_1_pz;
        MatX pd12_pz = ptc2_point_1_pz - ptc1_point_1_pz;

        // Solve for gradients of terms used to find distance
        VecX pa_pz = 2*batchDot(d1,pd1_pz);
        VecX pb_pz = batchDot(d1,pd2_pz)+batchDot(d2, pd1_pz);
        VecX pc_pz = batchDot(d12,pd1_pz)+batchDot(d1, pd12_pz);
        VecX pe_pz = 2*batchDot(d2,pd2_pz);
        VecX pf_pz = batchDot(d2,pd12_pz)+batchDot(d12, pd2_pz);
        VecX pg_pz = 2*batchDot(d12,pd12_pz);

        // auto time_pa = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time_pa - time1);
        // std::cout << "Time taken to find p_alphabet_pz: " << duration.count()/1e3 << " microseconds" << std::endl;

        // Solve for more complex terms
        VecX lgrey = pa_pz * tStar - pb_pz * uStar - pc_pz;
        auto dist_vec = d1 * tStar - d2 * uStar - d12;
        double purple = (dist_vec).norm();
        double orange = a * tStar - b * uStar - c;
        double green = -b * tStar - e * uStar + f;
        VecX red = -pb_pz * tStar + pe_pz * uStar + pf_pz;
        VecX grey = (batchDot((d1 * tStar*tStar - d2 * tStar * uStar - d12 * tStar),pd1_pz)-
                    batchDot((d1 * tStar * uStar - d2 * uStar*uStar - d12 * uStar),pd2_pz)-
                    batchDot((d1 - d2 - d12), pd12_pz))/ (purple*purple);

        // auto time_colors = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time_colors - time_pa);
        // std::cout << "Time taken to find colors: " << duration.count()/1e3 << " microseconds" << std::endl;

        // 
        double gradient_factor = 1/pow(purple,3);
        double L11;
        double L22;
        double L12;

        VecX dLambda1 = VecX::Zero(colmns);
        if (tStar == 1){
            dLambda1 = -(lgrey * purple - orange * grey);
        }
        VecX dLambda2 = VecX::Zero(colmns);
        if (tStar == 0){
            dLambda2 = -(red * purple - green * grey);
        }
        VecX dLambda3 = VecX::Zero(colmns);
        if (uStar == 1){
            dLambda3 = (lgrey * purple - orange * grey);
        }
        VecX dLambda4 = VecX::Zero(colmns);
        if (uStar == 0){
            dLambda4 = (red * purple - green * grey);
        }

        auto H_inner = batchDot(-d1 * tStar + d2 * uStar + d12,-pd1_pz * tStar + pd2_pz * uStar + pd12_pz);
        VecX H1;
        VecX H2;


        VecX pt_pz = VecX::Zero();
        VecX pu_pz = VecX::Zero();
        switch(ind){
            case 0:
            case 1:
            case 2:
            case 3:
                H2 = (
                    (red) / purple  + dLambda2 - dLambda4
                    - green * H_inner*gradient_factor
                ).transpose();
                L22 = e / purple - gradient_factor*green*green;
                pu_pz = -H2 / L22;
                break;
            case 4:
            case 5:
            case 6:
            case 7:
                H1 = (
                    (lgrey) / purple + dLambda1 - dLambda3
                    - orange * H_inner*gradient_factor
                ).transpose();
                L11 = a / purple - gradient_factor*orange*orange;
                pt_pz = -H1 / L11;
                break;
            case 8:
            case 9:
                H2 = (
                    (red) / purple  + dLambda2 - dLambda4
                    - green * H_inner*gradient_factor
                );
                H1 = (
                    (lgrey) / purple + dLambda1 - dLambda3
                    - orange * H_inner*gradient_factor
                );
                L11 = a / purple - gradient_factor*orange*orange;
                L22 = e / purple - gradient_factor*green*green;
                L12 = -b / purple - gradient_factor*orange * green;
                pt_pz = (L12 * H2 - L22 * H1) / (L11 * L22 - L12*L12);
                pu_pz = (L12 * H2 - L22 * H1) / (L11 * L22 - L12*L12);
                break;
            case 10:
            case 11:
            case 12:
            case 13:

                break;
        }

        // Calculate gradient of distance
        VecX part1 = batchDot(dist_vec, pd1_pz*tStar + d1*pt_pz.transpose() - pd2_pz*uStar - d2*pu_pz.transpose() - pd12_pz) / purple;

        pdist_pz = (part1 - r1 * pt_pz - r2 * pu_pz);

        // auto time2 = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time_colors);
        // std::cout << "Time taken to find grad: " << duration.count()/1e3 << " microseconds" << std::endl;
    }
    

    // auto stop = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    // std::cout << "Time taken by distInternal: " << duration.count()/1e3 << " microseconds" << std::endl;

    return finalDistance;
}

double TaperedCapsuleCollision::computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const MatX& ptc1_point_1_pz, const MatX& ptc1_point_2_pz, 
                        const MatX& ptc2_point_1_pz, const MatX& ptc2_point_2_pz, 
                        const double& tc1_radius_1, const double& tc1_radius_2, 
                        const double& tc2_radius_1, const double& tc2_radius_2,
                        VecX& pdist_pz){
    // auto start = std::chrono::high_resolution_clock::now();
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                ptc1_point_1_pz, ptc1_point_2_pz, ptc2_point_1_pz, ptc2_point_2_pz,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, pdist_pz);
                auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    // std::cout << "Time taken by computeDist: " << duration.count()/1e3 << " microseconds" << std::endl;
    return distance;
}

double TaperedCapsuleCollision::computeDistance(const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
                        const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
                        const double& tc1_radius_1, const double& tc1_radius_2, 
                        const double& tc2_radius_1, const double& tc2_radius_2){
    MatX nullDerivative;
    VecX nullGradient;
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                nullDerivative, nullDerivative, nullDerivative, nullDerivative,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, nullGradient);
    return distance;
}


} // namespace RAPTOR