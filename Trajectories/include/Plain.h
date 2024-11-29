
#ifndef PLAIN_H
#define PLAIN_H

#include "Trajectories.h"

namespace RAPTOR {

class Plain : public Trajectories {
public:
    // Constructor
    Plain() = default;

    // Constructor
    Plain(const int Nact_input);

    // Constructor
    Plain(const int N_input, const int Nact_input);

    // Destructor
    ~Plain() = default;
    
    // class methods:
        // compute trajectories
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override;
};

}; // namespace RAPTOR

#endif // PLAIN_H
