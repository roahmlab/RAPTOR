#include "TaperedCapsuleCollision.h"
#include <chrono>

using namespace RAPTOR;

int main() {
    // macro NUM_THREADS should be define in cmake
    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        throw std::runtime_error("macro NUM_THREADS is not defined!");
    #endif

    std::srand(time(NULL));
    const size_t NUM_GRAD = 7;
    const size_t instances = 64 * 9;
    double dist[instances] = {0};
    double total_time[instances] = {0};

    TaperedCapsuleCollision<NUM_GRAD> collider;

    Eigen::Vector3d p11 = Eigen::Vector3d::Random() * 5;
    Eigen::Vector3d p12 = Eigen::Vector3d::Random() * 5;
    Eigen::Vector3d p21 = Eigen::Vector3d::Random() * 5;
    Eigen::Vector3d p22 = Eigen::Vector3d::Random() * 5;

    Eigen::Matrix<double, 3, NUM_GRAD> p11_pz;
    Eigen::Matrix<double, 3, NUM_GRAD> p12_pz;
    Eigen::Matrix<double, 3, NUM_GRAD> p21_pz;
    Eigen::Matrix<double, 3, NUM_GRAD> p22_pz;
    p11_pz.setRandom();
    p12_pz.setRandom();
    p21_pz.setRandom();
    p22_pz.setRandom();
    
    double r11(0.5);
    double r12(1.0);
    double r21(0.5);
    double r22(1.0);
    
    Eigen::Vector<double, NUM_GRAD> dist_grad;
    auto begin = std::chrono::high_resolution_clock::now();
    int i = 0;
    #pragma omp parallel for shared(collider, p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad, dist) private(i) schedule(dynamic)
    for(i = 0; i < instances; i++){
        dist[i] = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - begin);
    std::cout << "Time taken total in parallel: " << duration.count() / 1.0e3 << " microseconds" << std::endl;

    dist[0] = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
    begin = std::chrono::high_resolution_clock::now();
    for(i = 0; i < instances; i++){
        dist[i] = collider.computeDistance(p11, p12, p21, p22, p11_pz, p12_pz, p21_pz, p22_pz, r11, r12, r21, r22, dist_grad);
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - begin);
    std::cout << "Time taken on average: " << duration.count() / (instances * 1.0e3) << " microseconds" << std::endl;

    return 0;
}