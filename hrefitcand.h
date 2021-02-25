#ifndef HREFITCAND_H
#define HREFITCAND_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <hvirtualcand.h>


class HRefitCand : public TLorentzVector {
private:
    Double_t fR, fZ; //fTheta, fPhi;
    TMatrixD fCov;
    HVirtualCand * cand;

public:
    HRefitCand(HVirtualCand * cand);
    void setR(Double_t val) { fR = val; }
    void setZ(Double_t val) { fZ = val; }
    //void setTheta(Double_t val) { fTheta = val; }
    //void setPhi(Double_t val) { fPhi = val; }

    void setCovariance(const TMatrixD & cov);
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    //Double_t getTheta() const { return fTheta; }
    //Double_t getPhi() const { return fPhi; }
    TMatrixD getCovariance() const {return fCov; }

    void reset();
    void update();
};

#endif /* HREFITCAND_H */
