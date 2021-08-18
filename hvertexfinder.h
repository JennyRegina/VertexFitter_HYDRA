/**
 * HVertexFinder.h
 *
 *
 */

#ifndef HVERTEXFINDER_H
#define HVERTEXFINDER_H

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

class HVertexFinder

 
{
private:
    std::vector<HRefitCand> fCands;

    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    int fVerbose;

    double fDistanceParticleToParticle;
    double fDistanceParticleToVertex;

    double fDistParticle1Vertex;
    double fDistParticle2Vertex;
    double fDistParticle1Origin;
    double fDistParticle2Origin;

    bool fPrimaryVertexFound;

    // Properties related to difference between primary and decay vertex
    TVector3 fVecPrimToDecayVertex;
    double fDistPrimToDecayVertex;

    // Checks if the decay vertex is located further downstream in z-position
    // and at larger values of R compared to the primary vertex
    bool fPrimaryVertexIsBetforeDecayVertex;
    bool fPrimaryVertexIsInsideDecayVertex;

    bool fUsePrimaryVertexInNeutralCandidateCalculation;

public:
    HVertexFinder();
    ~HVertexFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    TVector3 findVertex(const std::vector<HRefitCand> &cands);
    TVector3 findPrimaryVertex(const std::vector<HRefitCand> &cands);
    std::vector<HRefitCand> UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &VertexPos);
    void calculateVertexProperties(TVector3 primaryVertex, TVector3 decayVertex);

    TVector3 getVertex() const { return fVertex; }               // Function that the user should use in the analysis macro
    TVector3 gePrimarytVertex() const { return fPrimaryVertex; } // Function that the user should use in the analysis macro

    double getDistanceBetweenFittedParticles() const { return fDistanceParticleToParticle; }
    double getDistanceFirstParticleVertex() const { return fDistParticle1Vertex; }
    double getDistanceSecondParticleVertex() const { return fDistParticle2Vertex; }
    double getDistanceFirstParticleOrigin() const { return fDistParticle1Origin; }
    double getDistanceSecondParticleOrigin() const { return fDistParticle2Origin; }

    double getDistBetweenVertices() { return fDistPrimToDecayVertex; }

    bool isPrimVertexBeforeDecayVertex() { return fPrimaryVertexIsBetforeDecayVertex; }
    bool isPrimVertexInsideDecayVertex() { return fPrimaryVertexIsInsideDecayVertex; }

    TVector3 getPrimaryVertex() { return fPrimaryVertex; }
    
};

#endif /* HVERTEXFINDER_H */
