#include "ReachableSets.h"
#include "ArmourBezierCurves.h"
#include "InverseDynamics.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

int main() {
    // Test the PZsparse class
    PZsparse pz1(5.0);
    PZsparse pz2(3.0, 0.1f);
    
    // Test addition
    PZsparse pz3 = pz1 + pz2;
    
    // Test multiplication
    PZsparse pz4 = pz1 * 2.0;
    
    // Test division
    PZsparse pz5 = pz2 / 2.0;
    
    // Output results
    std::cout << "pz1: " << pz1 << std::endl;
    std::cout << "pz2: " << pz2 << std::endl;
    std::cout << "pz3 (pz1 + pz2): " << pz3 << std::endl;
    std::cout << "pz4 (pz1 * 2.0): " << pz4 << std::endl;
    std::cout << "pz5 (pz2 / 2.0): " << pz5 << std::endl;

    // Create a nontrivial PZsparse object
    // pz6 = 2 + 3 * k1 + 5 * qe1 ^ 2
    double pz6center = 2.0;
    double pz6coeff[] = {3.0, 5.0};
    uint32_t pz6degree[2][NUM_VARIABLES] = {0};
    pz6degree[0][0] = 1; // k1
    pz6degree[1][NUM_FACTORS] = 2; // qe1
    PZsparse pz6(pz6center, pz6coeff, pz6degree, 2);

    std::cout << pz6 << std::endl;

    // pz7 = pz6 * pz6
    //     = (2 + 3 * k1 + 5 * qe1 ^ 2) * (2 + 3 * k1 + 5 * qe1 ^ 2)
    //     = 4 + 12 * k1 + 20 * qe1 ^ 2 + 9 * k1 ^ 2 + 30 * k1 * qe1 ^ 2 + 25 * qe1 ^ 4
    double pz7center = 4.0;
    double pz7coeff[] = {12.0, 20.0, 9.0, 30.0, 25.0};
    uint32_t pz7degree[5][NUM_VARIABLES] = {0};
    pz7degree[0][0] = 1; // k1
    pz7degree[1][NUM_FACTORS] = 2; // qe1
    pz7degree[2][0] = 2; // k1 ^ 2
    pz7degree[3][0] = 1; // k1
    pz7degree[3][NUM_FACTORS] = 2; // k1 * qe1 ^ 2
    pz7degree[4][NUM_FACTORS] = 4; // qe1 ^ 4
    PZsparse pz7(pz7center, pz7coeff, pz7degree, 5);

    // verify pz7 = pz6 * pz6
    std::cout << pz7 - pz6 * pz6 << std::endl;

    // Create another nontrivial PZsparse object
    // pz8 = 2 + 0.01 * k1 ^ 2 + 0.02 * k4 + 0.05 * k2 ^ 2
    double pz8center = 2.0;
    double pz8coeff[] = {0.1, 0.2, 0.05};
    uint32_t pz8degree[3][NUM_VARIABLES] = {0};
    pz8degree[0][0] = 2; // k1
    pz8degree[1][3] = 1; // k4
    pz8degree[2][1] = 2; // k2
    PZsparse pz8(pz8center, pz8coeff, pz8degree, 3);

    // test sin and cos
    PZsparse pz9 = sin(pz8);
    PZsparse pz10 = cos(pz8);

    // verify sin^2 + cos^2 = 1
    PZsparse pz11 = pz9 * pz9 + pz10 * pz10;
    std::cout << pz11.toInterval() << std::endl;

    return 0;
}