// #include "InverseDynamics.h"
// #include "CustomizedInverseDynamics.h"
#include "TaperedCapsuleCollision.h"
#include <Eigen/Dense>

using namespace RAPTOR;

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

