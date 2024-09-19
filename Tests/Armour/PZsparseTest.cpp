#include "PZsparse.h"

using namespace RAPTOR;
using namespace Armour;

int main() {
    // Test the PZsparse class
    PZsparse pz1(5.0f);
    PZsparse pz2(3.0f, 0.1f);
    
    // Test addition
    PZsparse pz3 = pz1 + pz2;
    
    // Test multiplication
    PZsparse pz4 = pz1 * 2.0f;
    
    // Test division
    PZsparse pz5 = pz2 / 2.0f;
    
    // Output results
    std::cout << "pz1: " << pz1 << std::endl;
    std::cout << "pz2: " << pz2 << std::endl;
    std::cout << "pz3 (pz1 + pz2): " << pz3 << std::endl;
    std::cout << "pz4 (pz1 * 2.0): " << pz4 << std::endl;
    std::cout << "pz5 (pz2 / 2.0): " << pz5 << std::endl;

    return 0;
}