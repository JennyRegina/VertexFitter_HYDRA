// ROOT includes
#include "hrefitcand.h"

HRefitCand::HRefitCand(HVirtualCand* cand)
    : TLorentzVector(*cand), fR(cand->getR()), fZ(cand->getZ()), cand(cand)
{
}

void HRefitCand::setCovariance(const TMatrixD& cov)
{
    // 0 = 1/p
    // 1 = theta
    // 2 = phi
    // 3 = R
    // 4 = z
    fCov.ResizeTo(5, 5);
    fCov = cov;
}

void HRefitCand::reset() {
    *(TLorentzVector*)this = *cand;
}

void HRefitCand::update() {
    *((TLorentzVector*)cand) = *this;
}
