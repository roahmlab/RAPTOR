#define NUM_FACTORS 2

#include "TaperedCapsuleCollision.h"
#include <chrono>


using namespace RAPTOR;

int main() {
    std::srand(time(NULL));
    #define NUM_GRAD 7
    TaperedCapsuleCollision<NUM_GRAD> collider;
    int instances = 1000;
    double dist[instances] = {0};
    double total_time[instances] = {0};

    Eigen::Vector3d p11 = Eigen::Vector3d::Random()*5;
    Eigen::Vector3d p12 = Eigen::Vector3d::Random()*5;
    Eigen::Vector3d p21 = Eigen::Vector3d::Random()*5;
    Eigen::Vector3d p22 = Eigen::Vector3d::Random()*5;

    Eigen::Matrix<double,3,NUM_GRAD> p11_pz = Eigen::Matrix<double,3,7>::Random();
    Eigen::Matrix<double,3,NUM_GRAD> p12_pz = Eigen::Matrix<double,3,7>::Random();
    Eigen::Matrix<double,3,NUM_GRAD> p21_pz = Eigen::Matrix<double,3,7>::Random();
    Eigen::Matrix<double,3,NUM_GRAD> p22_pz = Eigen::Matrix<double,3,7>::Random();
    
    int i = 0;

    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);
    
    Eigen::Vector<double, NUM_GRAD> dist_grad;
    auto begin = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for num_threads(50)
    for(i = 0; i<instances; i++){
        
        auto start = std::chrono::high_resolution_clock::now();
        dist[i] = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        // std::cout << "Time taken by function: " << duration.count()/1e3 << " microseconds" << std::endl;
        // std::cout << "Avg Dist:" << dist.mean() << std::endl;
        total_time[i] = duration.count();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - begin);
    std::cout << "Time taken total: " << duration.count()/1.0e3 << " microseconds" << std::endl;
    double time_average = 0;
    for (int i = 0; i < instances; i++) {
        time_average += total_time[i];
    }
    time_average /= instances;
    std::cout << "Average function call: " << time_average/1.0e3 << " microseconds" << std::endl;

    return 0;
}