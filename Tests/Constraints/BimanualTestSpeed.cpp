#include "TaperedCapsuleCollision.h"
#include <chrono>

using namespace RAPTOR;

int main() {
    TaperedCapsuleCollision collider;
    Eigen::VectorXd dist(5000);
    auto start = std::chrono::high_resolution_clock::now();
    int i = 0;
    #pragma omp parallel for shared(dist) private(i) schedule(dynamic) num_threads(12)
    for(i = 0; i<5000; i++){
        Eigen::Vector3d p11 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p12 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p21 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p22 = Eigen::Vector3d::Random()*5;

        Eigen::MatrixXd p11_pz = Eigen::MatrixXd::Random(3,7);
        Eigen::MatrixXd p12_pz = Eigen::MatrixXd::Random(3,7);
        Eigen::MatrixXd p21_pz = Eigen::MatrixXd::Random(3,7);
        Eigen::MatrixXd p22_pz = Eigen::MatrixXd::Random(3,7);

        double r11(0.5);
        double r12(1.0);
        double r21(0.5);
        double r22(1.0);

        Eigen::MatrixXd dist_grad;
        dist(i) = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
        
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    std::cout << "Avg Dist:" << dist.mean() << std::endl;
    return 0;
}