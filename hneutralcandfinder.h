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

using std::cout;
using std::endl;

class HNeutralCandFinder
{
private:
    std::vector<HRefitCand> fCands;
    double fNeutralCandMass;

    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    int fVerbose=0;

    double fMomentumAfterDecay;

    HRefitCand fNeutralMotherCandidate;

    double fDistParticle1Vertex;
    double fDistParticle2Vertex;
    double fDistParticle1Origin;
    double fDistParticle2Origin;

    TMatrixD fCovarianceNeutralMother;
    bool fPrimaryVertexFound;

    double fPrimVtxResX;
    double fPrimVtxResY;
    double fPrimVtxResZ;

    double fDecVtxResX;
    double fDecVtxResY;
    double fDecVtxResZ;

public:
    HNeutralCandFinder(const std::vector<HRefitCand> &cands, double fNeutralCandMass);
    ~HNeutralCandFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // The first function is for creating a neutral mother candidate if only information of the decay vertex is available
    // The second function is for creating the neutral candidate if information about the primary vertex is also available
    void setNeutralMotherCand(double valMomentum, double valTheta, double valPhi, double valR, double ValZ, TVector3 decayVertex);
    void setNeutralMotherCand(TVector3 primVtx, TVector3 decayVtx);
    void setMassNutralCand(double val) { fNeutralCandMass = val; }

    void stetPrimaryVertexResolution(double valX, double valY, double valZ){ fPrimVtxResX=valX; fPrimVtxResY=valY; fPrimVtxResZ=valZ;}
    void stetDecayVertexResolution(double valX, double valY, double valZ){ fDecVtxResX=valX; fDecVtxResY=valY; fDecVtxResZ=valZ;}


    HRefitCand getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }
};

#endif /* HNEUTRALCANDFINDER_H */
