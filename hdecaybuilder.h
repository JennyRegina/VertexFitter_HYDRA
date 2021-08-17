/**
 * HDecayBuilder.h
 *
 *
 */

#ifndef HDECAYBUILDER_H
#define HDECAYBUILDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hcategorymanager.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"

#include "TH1F.h"

#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"

using std::cout;
using std::endl;

//double deg2rad = TMath::DegToRad();

class HDecayBuilder
{
private:
    // Working Particles
    std::vector<HRefitCand> fFitCands;
    std::vector< std::vector<HRefitCand> > fCands;
    // Output particles after fitting
    //std::vector<HFitParticleCand *> fOutputCands;
    std::vector<HRefitCand> fOutputCands;

    //Fitter input variables
    TString fTask;
    std::vector<Int_t> fPids;
    TLorentzVector fIniSys;
    HRefitCand fMother;
    Double_t fMass;
    
    Int_t fTotalCombos;
    Int_t fCombiCounter;
    std::vector<Int_t> particleCounter;
    
    Double_t fProb;

    int fVerbose;

public:
    //HDecayBuilder(std::vector<HRefitCand[]> fitCands, TString task, std::vector<Int_t> pids, TLorentzVector lv(0,0,0,0), HRefitCand mother=HRefitCand(), Double_t mass=0);
    HDecayBuilder(std::vector< std::vector<HRefitCand> > &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv = TLorentzVector(), HRefitCand mother = HRefitCand(), Double_t mass=0.);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }


    // Method to fill the data from a KParticleCand from old data
    //void FillData(KParticleCand * cand, HRefitCand & outcand, double arr[], double mass);
    
    //setters
    void setIniSys(TLorentzVector val) {fIniSys = val;}
    void setMother(HRefitCand val) {fMother = val;}
    void setMass(Double_t val) {fMass = val;}

    void buildDecay();
    
    void estimateCovarianceMatrix(HParticleCandSim *cand, HRefitCand *refitCand);

    void createNeutralCandidate();

    bool do4cFit();
    bool do3cFit();
    bool doMissMomFit();

    void fillFitCands();
    bool checkDoubleParticle(size_t i);

    // Functions for getting the pulls

    void createOutputParticle(HRefitCand FittedCand);
    void getFitCands(std::vector<HRefitCand> &cands) { cands = fOutputCands; }
    std::vector<HParticleCandSim> getOutput();

    void createOutputCategory();

    //void fillHistograms();
};

#endif /* HDECAYBUILDER_H */
