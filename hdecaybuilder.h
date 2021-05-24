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
    std::vector<HParticleCandSim *> fParticleCands;
    int fVerbose;

    // Variables used for setting the covariance matrix
    bool fFixedErrors;
    bool fMomDepErrors;

    // Containers for the protons, kaons and pions in the event
    std::vector<HRefitCand*> fProtons;
    std::vector<HRefitCand*> fPions;
    std::vector<HRefitCand*> fKaons;

    // Variables used for the vertex finding
    bool fFindPrimaryVertex;
    bool fFindDecayVertex;

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

    void addVertexFitting();

    void add3CFit();

};

#endif /* HDECAYBUILDER_H */
