#include "TaperedCapsuleCollision.h"
#include <chrono>

#define NUM_FACTORS 2

using namespace RAPTOR;

int main() {
    TaperedCapsuleCollision collider;
    int instances = 10000;
    Eigen::VectorXd dist(instances);
    double total_time = 0;
    
    int i = 0;
    for(i = 0; i<instances; i++){
        Eigen::Vector3d p11 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p12 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p21 = Eigen::Vector3d::Random()*5;
        Eigen::Vector3d p22 = Eigen::Vector3d::Random()*5;

        Eigen::Matrix<double,3,2> p11_pz = Eigen::Matrix<double,3,2>::Random();
        Eigen::Matrix<double,3,2> p12_pz = Eigen::Matrix<double,3,2>::Random();
        Eigen::Matrix<double,3,2> p21_pz = Eigen::Matrix<double,3,2>::Random();
        Eigen::Matrix<double,3,2> p22_pz = Eigen::Matrix<double,3,2>::Random();

        double r11(0.5);
        double r12(1.0);
        double r21(0.5);
        double r22(1.0);

        Eigen::Vector<double, 2> dist_grad;

        auto start = std::chrono::high_resolution_clock::now();
        dist(i) = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        // std::cout << "Time taken by function: " << duration.count()/1e3 << " microseconds" << std::endl;
        // std::cout << "Avg Dist:" << dist.mean() << std::endl;
        total_time = total_time*0.9+duration.count()*0.1;
    }
    std::cout << "Time taken by function: " << total_time/1e3 << " microseconds" << std::endl;

    return 0;
}