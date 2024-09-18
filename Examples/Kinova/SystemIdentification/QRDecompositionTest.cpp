// main.cpp

#include "QRDecompositionSolver.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>

int main() {

    int nq_nt = 100; 
    int p = 91;    

    RAPTOR::Kinova::QRDecompositionSolver qrSolver; 

    qrSolver.getData();


    std::cout << "Aid matrix:\n" << qrSolver.Aid << "\n\n";
    std::cout << "Ad matrix:\n" << qrSolver.Ad << "\n\n";
    std::cout << "Kd matrix:\n" << qrSolver.Kd << "\n\n";
    std::cout << "Beta vector:\n" << qrSolver.Beta << "\n\n";
    std::cout << "Ginv matrix:\n" << qrSolver.Ginv << "\n\n";
    std::cout << "Ginv RegroupMatrix:\n" << qrSolver.RegroupMatrix << "\n\n";
    std::cout << "Dimension of independent parameters: " << qrSolver.dim_id << "\n";
    std::cout << "Dimension of dependent parameters: " << qrSolver.dim_d << "\n";

  

    return 0;
}
