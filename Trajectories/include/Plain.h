
#ifndef PLAIN_H
#define PLAIN_H

#include "Trajectories.h"

namespace IDTO {

class Plain : public Trajectories {
public:
    // Constructor
    Plain() = default;

    // Constructor
    Plain(const int Nact_input);

    // Destructor
    ~Plain() = default;
    
    // class methods:
        // compute trajectories
    void compute(const VecX& z, bool compute_derivatives = true) override;
};

}; // namespace IDTO

#endif // PLAIN_H
