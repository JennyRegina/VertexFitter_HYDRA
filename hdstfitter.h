/**
 * HDSTFitter.h
 *
 *
 */

#ifndef HDSTFITTER_H
#define HDSTFITTER_H

// system includes
#include <iostream>
#include <vector>
#include <cmath>

// framework includes
#include "TLorentzVector.h"

#include "hdecaybuilder.h"
#include "hcovariancekinfit.h"

using std::cout;
using std::endl;

class HDSTFitter
{
private:
    TString fInfileList;
    bool fIncludeFw = false;
    std::vector< std::vector<HRefitCand> > fCandsFit;
    std:vector<Int_t> fPids;
    // Variables used for setting the covariance matrix
    bool fMomDepErrors = false;

    Int_t fVerbose = -1;

    // Method to fill the data from a HRefitCand for simulations
    void FillData(HParticleCandSim *cand, HRefitCand *outcand, double arr[], double mass);
    void FillDataFW(HFwDetCandSim *cand, HRefitCand *outcand, double arr[], double mass); //adjust to HForwardCand for newer Hydra
    //void estimateCov(Int_t pid, Double_t mom, double& cov[]); //double (&cov)[5]?

public:
    HDSTFitter(TString infilelist, bool includeFw, bool momDepErrors, Int_t nEvents);
    ~HDSTFitter(){};

    //User functions
    void addFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv = (0,0,0,0), Double_t mm = 0);
    void addBuilderTask(TString task, std::vector<Int_t> pids);

    void setIncludeFw(bool val){ fIncludeFw = val; }
    void setErrors();
    void setPids(std:vector<Int_t> val){ fPids = val; }
    void setVerbosity(Int_t val){ fVerbose = val; }

    std:vector<Int_t> getPids(){ return fPids; }


    void selectCandidates();
}

#endif /* HDSTFITTER_H */

// How to handle several tasks???
