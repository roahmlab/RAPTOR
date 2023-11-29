#include "KinovaCustomizedConstraints.h"

namespace IDTO {
namespace Kinova {

KinovaCustomizedConstraints::KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {
    m = 2 * trajPtr_->Nact;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void KinovaCustomizedConstraints::compute(const VecX& z, bool compute_derivatives) {
    trajPtr_->compute(z, compute_derivatives);
    
    g.head(trajPtr_->Nact) = trajPtr_->q(0);
    g.tail(trajPtr_->Nact) = trajPtr_->q(trajPtr_->N - 1);

    if (compute_derivatives) {
        pg_pz.topRows(trajPtr_->Nact) = trajPtr_->pq_pz(0);
        pg_pz.bottomRows(trajPtr_->Nact) = trajPtr_->pq_pz(trajPtr_->N - 1);
    }
}

void KinovaCustomizedConstraints::compute_bounds() {
    g_lb.head(trajPtr_->Nact).setConstant(-1);
    g_ub.head(trajPtr_->Nact).setConstant(-1);

    g_lb.tail(trajPtr_->Nact).setConstant(-1e19);
    g_ub.tail(trajPtr_->Nact).setConstant(1e19);
}

}; // namespace Kinova
}; // namespace IDTO