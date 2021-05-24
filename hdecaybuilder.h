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
    std::vector<HRefitCand> fCands;
    // Input particles
    std::vector<HParticleCandSim *> fParticleCands;
    // Output particles after fitting
    std::vector<HParticleCand *> fOutputCands;
    int fVerbose;

    // Variables used for setting the covariance matrix
    bool fFixedErrors=true;
    bool fMomDepErrors=false;

    // Containers for the protons, kaons and pions in the event
    std::vector<HRefitCand*> fProtons;
    std::vector<HRefitCand*> fPions;
    std::vector<HRefitCand*> fKaons;

    // Variables used for the vertex finding
    bool fFindPrimaryVertex=true;
    bool fFindDecayVertex=true;

public:
    HDecayBuilder(std::vector<HParticleCandSim*> particleCands);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // Method to fill the data from a HRefitCand for simulations
    void FillData(HParticleCandSim * cand, HRefitCand * outcand, double arr[], double mass);

    // Method to fill the data from a KParticleCand from old data
    //void FillData(KParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

    void estimateCovarianceMatrix(HParticleCandSim *cand, HRefitCand* refitCand);

    // Adds a vertex to the set of functions
    void addVertexFinding(bool valFindPrimaryVertex, bool valFindDecayVertex){

        fFindPrimaryVertex=valFindPrimaryVertex;
        fFindDecayVertex=valFindDecayVertex;

    };

    // Functions for the user to choose which fits to use
    void addVertexFit();
    void add3CFit();
    void add4CFit();
    void addMomentumFit();

    void createOutputParticle(HRefitCand);
    std::vector<HParticleCandSim> getOutput();

    void createOutputCategory();

};

#endif /* HDECAYBUILDER_H */
