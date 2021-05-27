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
    bool fFixedErrors = true;
    bool fMomDepErrors = false;

    bool fFitElectrons = false;

    // Containers for the particles in the event
    std::vector<HRefitCand *> fProtons;
    std::vector<HRefitCand *> fPions;
    std::vector<HRefitCand *> fKaons;
    std::vector<HRefitCand *> fNegKaons;
    std::vector<HRefitCand *> fPosPions;
    std::vector<HRefitCand *> fElectrons;
    std::vector<HRefitCand *> fPositrons;

    std::vector<HRefitCand> fCandsMissPos;

    // Variables used for the vertex finding
    bool fFindPrimaryVertex = true;
    bool fFindDecayVertex = true;

    // Variables used for the fitting
    bool fFitPrimaryVertex = true;
    bool fFitDecayVertex = true;
    bool fFit3C = true;
    bool fFit4C = false;
    bool fFitMomentum = false;

    // Probability cut values
    double fPrimVertexProbabilityCut = 0; // pow(10, -15);
    double fDecayVertexProbabilityCut = 0; //pow(10, -15);
    double f3CProbabilityCut = 0; // pow(10, -15);
    double f4CProbabilityCut = 0; // pow(10, -15);
    double fMomentumFitProbabilityCut = 0; // pow(10, -15);

    // Convergence values
    double fConvergenceCritPrimVtxFit = 0.01;
    double fConvergenceCritDecayVtxFit = 0.01;
    double fConvergenceCrit3CFit = 0.01;
    double fConvergenceCrit4CFit = 0.01;
    double fConvergenceCritMomentumFit = 0.01;

    // Maximum number of iterations
    int fIterationsPrimVtxFit = 10;
    int fIterationsDecayVtxFit = 10;
    int fIterations3CFit = 10;
    int fIterations4CFit = 10;
    int fIterationsMomentumFit = 10;

    // Chi2 values after the fit
    double fPrimVertexChi2 = -1;
    double fDecayVertexChi2 = -1;
    double f3CChi2 = -1;
    double f4CChi2 = -1;
    double fMomentumFitChi2 = -1;

    // Probability values as obtained from the fit
    double fPrimVertexProbability = -1;
    double fDecayVertexProbability = -1;
    double f3CProbability = -1;
    double f4CProbability = -1;
    double fMomentumFitProbability = -1;

    // Number of iterations of the different fits
    double fPrimVertexNumIter = -1;
    double fDecayVertexNumIter = -1;
    double f3CNumIter = -1;
    double f4CNumIter = -1;
    double fMomentumFitNumIter = -1;

    // Variables stating if the fit has converged
    bool fPrimVertexFitIsConverged = false;
    bool fDecayVertexFitIsConverged = false;
    bool f3CFitIsConverged = false;
    bool f4CFitIsConverged = false;
    bool fMomentumFitIsConverged = false;

public:
    HDecayBuilder(std::vector<HParticleCandSim *> particleCands);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // Method to fill the data from a HRefitCand for simulations
    void FillData(HParticleCandSim *cand, HRefitCand *outcand, double arr[], double mass);

    // Method to fill the data from a KParticleCand from old data
    //void FillData(KParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

    void estimateCovarianceMatrix(HParticleCandSim *cand, HRefitCand *refitCand);

    void createNeutralCandidate();

    void fitElectrons();

    void setFitElectrons(bool val)
    {
        fFitElectrons = val;
        if (val == true)
        {
            fMomDepErrors = true;
        }
    }

    // Adds a vertex to the set of functions
    void addVertexFinding(bool valFindPrimaryVertex, bool valFindDecayVertex)
    {
        fFindPrimaryVertex = valFindPrimaryVertex;
        fFindDecayVertex = valFindDecayVertex;
    };

    // Functions for the user to choose which fits to use
    void addVertexFit(bool val1, bool val2)
    {
        fFitPrimaryVertex = val1;
        fFitDecayVertex = val2;
    }
    void add3CFit(bool val) { fFit3C = val; }
    void add4CFit(bool val) { fFit4C = val; }
    void addMomentumFit(bool val) { fFitMomentum = val; }

    // User settings for the probability cuts
    void setPrimaryVertexProbabilityCut(double val) { fPrimVertexProbabilityCut = val; }
    void setDecayVertexProbabilityCut(double val) { fDecayVertexProbabilityCut = val; }
    void set3CProbabilityCut(double val) { f3CProbabilityCut = val; }
    void set4CProbabilityCut(double val) { f4CProbabilityCut = val; }
    void setMomentumFitProbabilityCut(double val) { fMomentumFitProbabilityCut = val; }

    // User settings for the convergence criteria of the different fits
    void setConvergenceCritPrimVtxFit(double val) { fConvergenceCritPrimVtxFit = val; }
    void setConvergenceCritDecayVtxFit(double val) { fConvergenceCritDecayVtxFit = val; }
    void setConvergenceCrit3CFit(double val) { fConvergenceCrit3CFit = val; }
    void setConvergenceCrit4CFit(double val) { fConvergenceCrit4CFit = val; }
    void setConvergenceCritMomentumFit(double val) { fConvergenceCritMomentumFit = val; }

    // User settings for the maximum number of iterations
    void setMaxIterationsPrimVtxFit(double val) { fIterationsPrimVtxFit = val; }
    void setMaxIterationsDecayVtxFit(double val) { fIterationsDecayVtxFit = val; }
    void setMaxIterations3CFit(double val) { fIterations3CFit = val; }
    void setMaxIterations4CFit(double val) { fIterations4CFit = val; }
    void setMaxIterationsMomentumFit(double val) { fIterationsMomentumFit = val; }

    // Functions for getting the chi2 of the fit
    double getPrimaryVertexChi2() { return fPrimVertexChi2; }
    double getDecayVertexChi2() { return fDecayVertexChi2; }
    double get3CChi2() { return f3CChi2; }
    double get4CChi2() { return f4CChi2; }
    double getMomentumChi2() { return fMomentumFitChi2; }

    // Functions for getting the probabilities of the fits
    double getPrimaryVertexProbability() { return fPrimVertexProbability; }
    double getDecayVertexProbability() { return fDecayVertexProbability; }
    double get3CProbability() { return f3CProbability; }
    double get4CProbability() { return f4CProbability; }
    double getMomentumFitProbability() { return fMomentumFitProbability; }

    // Functions for getting the number of iterations
    double getPrimaryVertexNumIter() { return fPrimVertexNumIter; }
    double getDecayVertexNumIter() { return fDecayVertexNumIter; }
    double get3CNumIter() { return f3CNumIter; }
    double get4CNumIter() { return f4CNumIter; }
    double getMomentumFitNumIter() { return fMomentumFitNumIter; }

    // Functions for obtaining information about if an event has converged
    bool primaryVertexIsConverged() { return fPrimVertexFitIsConverged; }
    bool decayVertexIsConverged() { return fDecayVertexFitIsConverged; }
    bool threeCFitIsConverged() { return f3CFitIsConverged; }
    bool fourCFitIsConverged() { return f4CFitIsConverged; }
    bool momentumFitIsConverged() { return fMomentumFitIsConverged; }

    // Functions for getting the pulls

    void createOutputParticle(HRefitCand);
    std::vector<HParticleCandSim> getOutput();

    void createOutputCategory();
};

#endif /* HDECAYBUILDER_H */