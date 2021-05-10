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
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"

using std::cout;
using std::endl;

class HDecayBuilder
{
private:
    std::vector<HRefitCand> fCands;
    int fVerbose;

    // Variables used for setting the covariance matrix
    bool fFixedErrors;

    // Containers for the protons, kaons and pions in the event
    std::vector<HRefitCand> fProtons;
    std::vector<HRefitCand> fPions;
    std::vector<HRefitCand> fKaons;

    // Variables used for the vertex finding
    bool fFindPrimaryVertex;
    bool fFindDecayVertex;

    // Variables used for setting the covariance matrix
    bool fFixedErrors;
    bool fMomDepErrors;

    std::vector<HRefitCand> 

public:
    HDecayBuilder();
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // Method to fill the data from a HRefitCand for simulations
    void FillData(HParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

    // Method to fill the data from a KParticleCand from old data
    void FillData(KParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

    // Adds a vertex to the set of functions
    void addVertexFinding(valFindPrimaryVertex,valFindDecayVertex){

        fFindPrimaryVertex=valFindPrimaryVertex;
        fFindDecayVertex=valFindDecayVertex;

    };

    void addVertexFitting();

    void add3CFit();

};

#endif /* HDECAYBUILDER_H */
