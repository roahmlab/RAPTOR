#include "TaperedCapsuleCollision.h"
#include <chrono>
#include <iomanip>

namespace RAPTOR {

inline void solve_quadratic(const double a, 
                            const double b, 
                            const double c,
                            double* sol1,
                            double* sol2) {
    double discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) {
        *sol1 = -1.0;
        *sol2 = -1.0;
        return;
    }
    double sqrt_discriminant = std::sqrt(discriminant);
    double twoa = 2 * a;
    *sol1 = (-b + sqrt_discriminant) / twoa;
    *sol2 = (-b - sqrt_discriminant) / twoa;
}

// explicit instantiation for cases currently used
template class TaperedCapsuleCollision<2>; // used in tests
template class TaperedCapsuleCollision<7>; // used in Armour
template class TaperedCapsuleCollision<14>; // used in DualArmour
template class TaperedCapsuleCollision<7 * 1 * 3 * 2>; // used in DualKinovaExample

template<int factors>
inline Eigen::Vector<double,factors> TaperedCapsuleCollision<factors>::batchDot(
    const Vec3 vector, 
    const Mat3F matrix) {
    // equivalent to this line, but for loop is faster
    // (matrix.transpose()*vector)
    Eigen::VectorXd result(factors);
    for(uint i = 0; i < factors; i++) {
        result[i] = vector.dot(matrix.col(i));
    }
    return result;
}

template<int factors> inline double TaperedCapsuleCollision<factors>::distanceInternal(
    const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
    const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
    const Mat3F& ptc1_point_1_pz, const Mat3F& ptc1_point_2_pz, 
    const Mat3F& ptc2_point_1_pz, const Mat3F& ptc2_point_2_pz, 
    const double tc1_radius_1, const double tc1_radius_2, 
    const double tc2_radius_1, const double tc2_radius_2,
    VecF& pdist_pz) {

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

    double distances[14];
    double u_test[14];
    double t_test[14];
    
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
    if(std::abs(r2) > 1e-10){
        solve_quadratic(e*(e-r2*r2), 2*f*(e-r2*r2), f*f-r2*r2*g, &u_test[0], &u_test[1]);
        solve_quadratic(e*(e-r2*r2), 2*(f-b)*(e-r2*r2), (f-b)*(f-b)-r2*r2*(g-2*c+a), &u_test[2], &u_test[3]);    
    }
    else{
        u_test[0] = -f*(e-r2*r2)/(e*(e-r2*r2));
        u_test[1] = -1;
        u_test[2] = -(f-b)*(e-r2*r2)/(e*(e-r2*r2));
        u_test[3] = -1;
    }

    t_test[0] = 0;
    t_test[1] = 0;
    t_test[2] = 1;
    t_test[3] = 1;

    u_test[4] = 0;
    u_test[5] = 0;
    u_test[6] = 1;
    u_test[7] = 1;

        // case 3-4
    if(std::abs(r1) > 1e-10){
        solve_quadratic(a*(a-pow(r1,2)), -2*c*(a-pow(r1,2)), pow(c,2)-pow(r1,2)*g, &t_test[4], &t_test[5]);
        solve_quadratic(a*(a-pow(r1,2)), -2*(b+c)*(a-pow(r1,2)), pow((b+c),2)-pow(r1,2)*(e+2*f+g), &t_test[6], &t_test[7]);
    }
    else{
        t_test[4] = c*(a-pow(r1,2))/(a*(a-pow(r1,2)));
        t_test[5] = -1;
        t_test[6] = (b+c)*(a-pow(r1,2))/(a*(a-pow(r1,2)));
        t_test[7] = -1;
    }

    
    
    // edge to edge check
        // case 5
    double radicand = (g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e);
    if(radicand < 0) {
        t_test[8] = -1.0;
        t_test[9] = -1.0;
        u_test[8] = -1.0;
        u_test[9] = -1.0;
    }
    else{
        if(radicand < 0) {
            radicand = 0;
        }
        double sqrt_radicand = sqrt(radicand);
        t_test[8] = (c*r2 + f*r1)/(a*r2 + b*r1) + ((b*r2 + e*r1)*(a*r2*sqrt_radicand - b*b*b*c + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));
        t_test[9] = (c*r2 + f*r1)/(a*r2 + b*r1) - ((b*r2 + e*r1)*(b*b*b*c + a*r2*sqrt_radicand + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
         - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));
        u_test[8] = (a*r2*sqrt_radicand - b*b*b*c + b*r1*sqrt_radicand
         + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);
        u_test[9] = -(b*b*b*c + a*r2*sqrt_radicand + b*r1*sqrt_radicand
         - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);
    }    

    for(int i = 0; i < 14; i++) {
        if(u_test[i] <= 1 && u_test[i] >= 0 && 
           t_test[i] <= 1 && t_test[i] >= 0) {
            distances[i] = (d1*t_test[i]-d2*u_test[i]-d12).norm() - r1*t_test[i] - r2*u_test[i];
        }
        else{
            distances[i] = INFINITY;
        }
    }

    size_t ind = 0;
    for (int i = 1; i < 14; i++) {
        if (distances[i] < distances[ind]) {
            ind = i;
        }
    }

    double finalDistance = distances[ind] - tc1_radius_1 - tc2_radius_1;
    double uStar = u_test[ind];
    double tStar = t_test[ind];

    if(ptc1_point_1_pz.size()>0) {
        // Solve for gradients of T and U
        Mat3F pd1_pz = ptc1_point_2_pz - ptc1_point_1_pz;
        Mat3F pd2_pz = ptc2_point_2_pz - ptc2_point_1_pz;
        Mat3F pd12_pz = ptc2_point_1_pz - ptc1_point_1_pz;

        // Solve for gradients of terms used to find distance
        VecF pa_pz = 2 * batchDot(d1, pd1_pz);
        VecF pb_pz = batchDot(d1, pd2_pz) + batchDot(d2, pd1_pz);
        VecF pc_pz = batchDot(d12, pd1_pz) + batchDot(d1, pd12_pz);
        VecF pe_pz = 2 * batchDot(d2, pd2_pz);
        VecF pf_pz = batchDot(d2, pd12_pz) + batchDot(d12, pd2_pz);
        VecF pg_pz = 2 * batchDot(d12, pd12_pz);

        // Solve for more complex terms
        VecF lgrey = pa_pz * tStar - pb_pz * uStar - pc_pz;
        auto dist_vec = d1 * tStar - d2 * uStar - d12;
        double purple = dist_vec.norm();
        // To reach zero the capsules are already in collision, but this protects against undefined gradient
        if(purple==0) {
            purple = 1e-6;
        }
        double orange = a * tStar - b * uStar - c;
        double green = -b * tStar - e * uStar + f;
        VecF red = -pb_pz * tStar + pe_pz * uStar + pf_pz;
        auto grey_intermediate = d1 * tStar * tStar - d2 * tStar * uStar - d12 * tStar;
        VecF grey = (batchDot(grey_intermediate, pd1_pz) -
                     batchDot(grey_intermediate, pd2_pz) -
                     batchDot((d1 - d2 - d12), pd12_pz)) / (purple*purple);

        double gradient_factor = 1 / pow(purple, 3);
        double L11;
        double L22;
        double L12;

        VecF dLambda1 = VecF::Zero();
        if (tStar == 1) {
            dLambda1 = -(lgrey * purple - orange * grey);
        }
        VecF dLambda2 = VecF::Zero();
        if (tStar == 0) {
            dLambda2 = -(red * purple - green * grey);
        }
        VecF dLambda3 = VecF::Zero();
        if (uStar == 1) {
            dLambda3 = (lgrey * purple - orange * grey);
        }
        VecF dLambda4 = VecF::Zero();
        if (uStar == 0) {
            dLambda4 = (red * purple - green * grey);
        }

        auto H_inner = batchDot(-d1 * tStar + d2 * uStar + d12,-pd1_pz * tStar + pd2_pz * uStar + pd12_pz);

        VecF H1;
        VecF H2;
        VecF pt_pz = VecF::Zero();
        VecF pu_pz = VecF::Zero();;

        switch(ind) {
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
        VecF part1 = batchDot(dist_vec, pd1_pz * tStar + d1 * pt_pz.transpose() - pd2_pz * uStar - d2 * pu_pz.transpose() - pd12_pz) / purple;

        pdist_pz = (part1 - r1 * pt_pz - r2 * pu_pz);
    }
    
    return finalDistance;
}

template<int factors> double TaperedCapsuleCollision<factors>::computeDistance(
    const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
    const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
    const Mat3F& ptc1_point_1_pz, const Mat3F& ptc1_point_2_pz, 
    const Mat3F& ptc2_point_1_pz, const Mat3F& ptc2_point_2_pz, 
    const double tc1_radius_1, const double tc1_radius_2, 
    const double tc2_radius_1, const double tc2_radius_2,
    VecF& pdist_pz) {
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                ptc1_point_1_pz, ptc1_point_2_pz, ptc2_point_1_pz, ptc2_point_2_pz,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, pdist_pz);
                auto stop = std::chrono::high_resolution_clock::now();
    return distance;
}

template<int factors> double TaperedCapsuleCollision<factors>::computeDistance(
    const Vec3& tc1_point_1, const Vec3& tc1_point_2, 
    const Vec3& tc2_point_1, const Vec3& tc2_point_2, 
    const double tc1_radius_1, const double tc1_radius_2, 
    const double tc2_radius_1, const double tc2_radius_2) {
    Mat3F nullDerivative;
    VecF nullGradient;
    double distance = distanceInternal(tc1_point_1, tc1_point_2, tc2_point_1, tc2_point_2, 
                nullDerivative, nullDerivative, nullDerivative, nullDerivative,
                tc1_radius_1, tc1_radius_2, tc2_radius_1, tc2_radius_2, nullGradient);
    return distance;
}


} // namespace RAPTOR