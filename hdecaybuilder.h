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

#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"

using std::cout;
using std::endl;

double deg2rad = TMath::DegToRad();

class HDecayBuilder
{
private:
    // Working Particles
    std::vector<HRefitCand[]> fCands;
    // Output particles after fitting
    std::vector<HFitParticleCand *> fOutputCands;
    int fVerbose;

    //Fitter input variables
    std::vector<Int_t> fPids;
    TLorentzVector fIniSys;
    HRefitCand fMother;
    Double_t fMass;


public:
    HDecayBuilder(std::vector<HRefitCand[]> fitCands, TString task, std::vector<Int_t> pids, TLorentzVector lv = (0,0,0,0), HRefitCand mother, Double_t mass=0);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }


    // Method to fill the data from a KParticleCand from old data
    //void FillData(KParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

    void estimateCovarianceMatrix(HParticleCandSim *cand, HRefitCand *refitCand);

    void createNeutralCandidate();

    void do4cFit();
    void do3cFit();
    void doMissMomentumFit();

    void fillFitCands();

    void buildDecay();

    // Functions for getting the pulls

    void createOutputParticle(HRefitCand);
    std::vector<HParticleCandSim> getOutput();

    void createOutputCategory();
};

#endif /* HDECAYBUILDER_H */