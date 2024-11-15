// #include "InverseDynamics.h"
// #include "CustomizedInverseDynamics.h"
#include "TaperedCapsuleCollision.h"

using namespace RAPTOR;

// float solve_quadratic(float a, float b, float c, int sign){
//     return (-b+sign*sqrt(pow(b,2)-4*a*c))/(2*a);
// }

// float computeDistance(const taperedCapsule tc1, const taperedCapsule tc2){
//     using Vec3 = Eigen::Vector3f;
//     using VecX = Eigen::VectorXf;
//     using MatX = Eigen::MatrixXf;
//     std::cout << "test\n";
//     float r1 = tc1.r1-tc1.r2;
//     float r2 = tc2.r1-tc2.r2;

//     Vec3 d1 = tc1.p2-tc1.p1;
//     Vec3 d2 = tc2.p2-tc2.p1;
//     Vec3 d12 = tc2.p1-tc1.p1;

//     std::cout << "D1\n" << d1 << "\n";

//     float a = d1.dot(d1);
//     float b = d1.dot(d2);
//     float c = d1.dot(d12);
//     float e = d2.dot(d2);
//     float f = d2.dot(d12);
//     float g = d12.dot(d12);

//     std::cout << "a\n" << a << "\n";

//     Eigen::VectorXf distances(14);
//     Eigen::VectorXf u_test(14);
//     Eigen::VectorXf t_test(14);
    
//     // end point to end point checks
//     u_test[10] = 0.0;
//     u_test[11] = 0.0;
//     u_test[12] = 1.0;
//     u_test[13] = 1.0;
//     t_test[10] = 0.0;
//     t_test[11] = 1.0;
//     t_test[12] = 0.0;
//     t_test[13] = 1.0;

//     std::cout << "u_test\n" << u_test[0] << "\n";
    
//     // end point to edge checks
//         // case 1-2
//     u_test[0] = solve_quadratic(e*(e-r2*r2), 2*f*(e-r2*r2), f*f-r2*r2*g, 1);
//     u_test[1] = solve_quadratic(e*(e-r2*r2), 2*f*(e-r2*r2), f*f-r2*r2*g, -1);
//     u_test[2] = solve_quadratic(e*(e-r2*r2), 2*(f-b)*(e-r2*r2), (f-b)*(f-b)-r2*r2*(g-2*c+a), 1);
//     u_test[3] = solve_quadratic(e*(e-r2*r2), 2*(f-b)*(e-r2*r2), (f-b)*(f-b)-r2*r2*(g-2*c+a), -1);
//     t_test[0] = 0;
//     t_test[1] = 0;
//     t_test[2] = 1;
//     t_test[3] = 1;
//         // case 3-4
//     u_test[4] = 0;
//     u_test[5] = 0;
//     u_test[6] = 1;
//     u_test[7] = 1;
//     t_test[4]  = solve_quadratic(a*(a-pow(r1,2)), -2*c*(a-pow(r1,2)), pow(c,2)-pow(r1,2)*g, 1);
//     t_test[5]  = solve_quadratic(a*(a-pow(r1,2)), -2*c*(a-pow(r1,2)), pow(c,2)-pow(r1,2)*g, -1);
//     t_test[6] = solve_quadratic(a*(a-pow(r1,2)), -2*(b+c)*(a-pow(r1,2)), pow((b+c),2)-pow(r1,2)*(e+2*f+g), 1);
//     t_test[7] = solve_quadratic(a*(a-pow(r1,2)), -2*(b+c)*(a-pow(r1,2)), pow((b+c),2)-pow(r1,2)*(e+2*f+g), -1);
//     // edge to edge check
//         // case 5
//     float phi = (r2*b+r1*f)/(r2*a+r1*b);
//     float gamma = (r2*c+r1*f)/(r2*a+r1*b);
//     float alpha = a*phi-b;
//     float beta = a*gamma-c;
    
//     float a_5 = alpha*alpha - r1*r1*phi*phi*a - r1*r1*e + 2*r1*r1*b*phi;
//     float b_5 = 2*alpha*beta - 2*r1*r1*phi*gamma*a + 2*r1*r1*b*gamma + 2*r1*c*phi - 2*r1*r1*f;
//     float c_5 = beta*beta - r1*r1*a*gamma*gamma + 2*r1*r1*c*gamma - r1*r1*g;
//     u_test[8] = solve_quadratic(a_5, b_5, c_5, 1);
//     u_test[9] = solve_quadratic(a_5, b_5, c_5, -1);
//     t_test[8] = phi*u_test[8]+gamma;
//     t_test[9] = phi*u_test[9]+gamma;
//     std::cout << u_test[8] << "\n" << u_test[9] << "\n";
//     std::cout << t_test[8] << "\n" << t_test[9] << "\n";
//     t_test[8] = (c*r2 + f*r1)/(a*r2 + b*r1) + ((b*r2 + e*r1)*(a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*b*b*c + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
//          + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));
//     t_test[9] = (c*r2 + f*r1)/(a*r2 + b*r1) - ((b*r2 + e*r1)*(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
//          - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2))/((a*r2 + b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1));
//     //t_test[10] = (c*r2 - f*r1)/(a*r2 - b*r1) - ((b*r2 - e*r1)*(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 - 2*b*b*c*r1*r2 - a*b*c*e + 2*a*b*f*r1*r2))/((a*r2 - b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 + 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b - 2*b*b*b*r1*r2 + b*b*e*r1*r1));
//     //t_test[11] = (c*r2 - f*r1)/(a*r2 - b*r1) - ((b*r2 - e*r1)*(b*b*b*c - a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 - 2*b*b*c*r1*r2 - a*b*c*e + 2*a*b*f*r1*r2))/((a*r2 - b*r1)*(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 + 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b - 2*b*b*b*r1*r2 + b*b*e*r1*r1));
    
//     u_test[8] = (a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*b*b*c + b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
//          + a*a*f*r2*r2 + a*b*b*f - a*a*e*f - a*b*c*r2*r2 - b*c*e*r1*r1 + a*e*f*r1*r1 - 2*b*b*c*r1*r2 + a*b*c*e + 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);
//     u_test[9] = -(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b + 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e))
//          - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 + 2*b*b*c*r1*r2 - a*b*c*e - 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 - 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b + 2*b*b*b*r1*r2 + b*b*e*r1*r1);
//     std::cout << "matlab\n" << u_test[8] << "\n" << u_test[9] << "\n";
//     std::cout << t_test[8] << "\n" << t_test[9] << "\n";
//     // u_test[10] = -(b*b*b*c + a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 - 2*b*b*c*r1*r2 - a*b*c*e + 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 + 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b - 2*b*b*b*r1*r2 + b*b*e*r1*r1);
//     // u_t  est[11] = -(b*b*b*c - a*r2*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) + b*r1*sqrt((g*b*b - 2*b*c*f + e*c*c + a*f*f - a*e*g)*(b*b - 2*b*r1*r2 + e*r1*r1 + a*r2*r2 - a*e)) - a*a*f*r2*r2 - a*b*b*f + a*a*e*f + a*b*c*r2*r2 + b*c*e*r1*r1 - a*e*f*r1*r1 - 2*b*b*c*r1*r2 - a*b*c*e + 2*a*b*f*r1*r2)/(a*a*e*e - a*a*e*r2*r2 - 2*a*b*b*e + a*b*b*r2*r2 + 2*a*b*e*r1*r2 - a*e*e*r1*r1 + b*b*b*b - 2*b*b*b*r1*r2 + b*b*e*r1*r1);
    
    

//     std::cout << "u_test\n" << u_test << "\n";
//     std::cout << "t_test\n" << t_test << "\n";

//     for(int i = 0; i<14; i++){
//         if(u_test[i] <= 1 && u_test[i] >= 0 && t_test[i] <= 1 && t_test[i] >= 0){
//             distances[i] = (d1*t_test[i]-d2*u_test[i]-d12).norm() - r1*t_test[i] - r2*u_test[i] - tc1.r1 - tc2.r1;
//         }
//         else{
//             distances[i] = INFINITY;
//         }
//     }
//     std::cout << "distances\n" << distances << "\n";
//     return distances.minCoeff();
// }

int main() {
    Eigen::Vector3d p11(1.0, 0.0, 1.0);
    Eigen::Vector3d p12(0.0, 0.0, 0.0);
    Eigen::Vector3d p21(1.0, 3.0, -1.0);
    Eigen::Vector3d p22(0.0, 4.0, 0.0);

    double r11(1.0);
    double r12(1.5);
    double r21(1.0);
    double r22(1.5);

    double distance = 0.0;

    // std::unique_ptr<TaperedCapsuleCollision> contraintPtr;
    // contraintPtr = std::make_unique<TaperedCapsuleCollision>();
    TaperedCapsuleCollision collider;
    // std::cout << p11 << "\n";
    distance = collider.computeDistance(p11, p12, p21, p22, r11, r12, r21, r22);

    std::cout << "distance\n" << distance << "\n";
}

