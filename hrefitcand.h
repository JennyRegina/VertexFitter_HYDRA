#ifndef HREFITCAND_H
#define HREFITCAND_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include "/lustre/hades/user/wesmail/forJana/KParticleCand.h"


class HRefitCand : public TLorentzVector {
private:
    Double_t fR, fZ;
    TMatrixD fCov;
    KParticleCand * cand;

public:
  //  HRefitCand();
    HRefitCand(KParticleCand * cand);
    void setR(Double_t val) { fR = val; }
    void setZ(Double_t val) { fZ = val; }
    void setCovariance(const TMatrixD & cov);
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    TMatrixD getCovariance() const {return fCov; }

    void reset();
    void update();
};

#endif /* HREFITCAND_H */
