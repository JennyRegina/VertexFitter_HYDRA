#ifndef HREFITCAND_H
#define HREFITCAND_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

#include "/lustre/hades/user/jrieger/pp35_data_4charged/forJana/KParticleCand.h"

class HRefitCand : public TLorentzVector {
private:
    KParticleCand * cand;
    Double_t fR, fZ; //fTheta, fPhi;
    TMatrixD fCov;

public:
    HRefitCand(KParticleCand * cand);
    HRefitCand(); //Does that work like this? See https://os.mbed.com/users/fpucher/code/HIM0Board/wiki/Vererbung-in-C%2B%2B
    ~HRefitCand(){};
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
