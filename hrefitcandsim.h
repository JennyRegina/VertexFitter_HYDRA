#ifndef HREFITCANDSIM_H
#define HREFITCANDSIM_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <hvirtualcand.h>
#include <hrefitcand.h>


class HRefitCandSim : public HRefitCand {
private:
    Int_t fTrack=-2, fParentTrack=-2;
    Int_t fPID=-2, fParentID=-2;
    Float_t fGeantxVertex=0, fGeantyVertex=0, fGeantzVertex=0;
    Float_t fGeantTotMom=0;

public:
    HRefitCandSim();
    ~HRefitCandSim(){};
    void setTrackGeant(Int_t val) { fTrack = val; }
    void setParentTrackGeant(Int_t val) { fParentTrack = val; }
    void setPID(Int_t val) { fPID = val; }
    void setParentID(Int_t val) { fParentID = val; }
    void setGeantxVertex(Int_t val) { fGeantxVertex = val; }
    void setGeantyVertex(Int_t val) { fGeantyVertex = val; }
    void setGeantzVertex(Int_t val) { fGeantzVertex = val; }
    void setGeantTotMom(Int_t val) { fGeantTotMom = val; }

    Int_t getTrackGeant() const { return fTrack; }
    Int_t getParentTrackGeant() const { return fParentTrack; }
    Int_t getPID() const { return fPID; }
    Int_t getParentID() const { return fParentID; }
    Int_t getGeantxVertex() const { return fGeantxVertex; }
    Int_t getGeantyVertex() const { return fGeantyVertex; }
    Int_t getGeantzVertex() const { return fGeantzVertex; }
    Int_t getGeantTotMom() const { return fGeantTotMom; }

};

#endif /* HREFITCAND_H */
