/**
 * HVertexFinder.h
 *
 *
 */

#ifndef HNEUTRALCANDFINDER_H
#define HNEUTRALCANDFINDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "hvertexfinder.h"
#include "hgeomvector.h"
#include "hparticletool.h"

//#include "/lustre/hades/user/jrieger/pp35_data_4charged/forJana/KParticleCand.h"

using std::cout;
using std::endl;

class HNeutralCandFinder
{
private: 
    std::vector<HRefitCand> fCands;

    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    int fVerbose;

    double fMomentumAfterDecay;
    double fNeutralCandMass;
    
    KParticleCand fNeutralMotherCandidate;

    double fDistParticle1Vertex;
    double fDistParticle2Vertex;
    double fDistParticle1Origin;
    double fDistParticle2Origin;

    TMatrixD fCovarianceNeutralMother;
    bool fPrimaryVertexFound;

    bool fUsePrimaryVertexInNeutralCandidateCalculation;

public:
    HNeutralCandFinder(const std::vector<HRefitCand> &cands);
    ~HNeutralCandFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // The first function is for creating a neutral mother candidate if only information of the decay vertex is available
    // The second function is for creating the neutral candidate if information about the primary vertex is also available
    void setNeutralMotherCand(double valMomentum, double valTheta, double valPhi, double valR, double ValZ, TVector3 decayVertex);
    void setNeutralMotherCandFromPrimaryVtxInfo(TVector3 primVtx, TVector3 decayVtx);
    void setUsePrimaryVertexInNeutralMotherCalculation(bool val){ fUsePrimaryVertexInNeutralCandidateCalculation = val; }

    KParticleCand getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }

};

#endif /* HNEUTRALCANDFINDER_H */
