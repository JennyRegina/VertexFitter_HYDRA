/**
 * HVertexFitter.h
 *
 *
 */

#ifndef HVERTEXFITTER_H
#define HVERTEXFITTER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"

using std::cout;
using std::endl;

const double pi2 = TMath::PiOver2();

template <typename T>
void Print(T const &matrix)
{
    Int_t nrows = matrix.GetNrows();
    Int_t ncols = matrix.GetNcols();

    cout << "shape(" << nrows << "," << ncols << ")" << endl;

    for (Int_t i = 0; i < nrows; i++)
    {
        for (Int_t j = 0; j < ncols; j++)
        {
            Double_t element = matrix(i, j);
            if (TMath::Abs(element) < 1e-10)
                element = 0.;
            if (element >= 0.)
                cout << " " << std::fixed << std::setw(8) << std::scientific << element << " ";
            else
                cout << std::fixed << std::setw(8) << std::scientific << element << " ";
        }
        cout << endl;
    }

    cout << endl;
}

class HVertexFitter
{
private:
    TMatrixD y, V, fPull;
    double fChi2, fProb;
    bool fConverged;
    int fIteration, fN;
    std::vector<HRefitCand> fCands;

    // data members for constraints
    double fMass;
    int fNdf;

    std::vector<double> fM;
    TLorentzVector fInit;
    bool fVtxConstraint;
    double fVtxPos;
    TVector3 fVertex;
    int fVerbose;

    double fLearningRate;
    int fNumIterations;

    double fDistanceParticleToParticle;
    double fDistanceParticleToVertex;

    double fPhi1Original;
    double fPhi2Original;
    HVirtualCand fLambdaCandidate;
public:
    HVertexFitter(const std::vector<HRefitCand> &cands);
    ~HVertexFitter(){};
    TMatrixD f_eval(const TMatrixD &m_iter);
    TMatrixD Feta_eval(const TMatrixD &miter);
    TVector3 findVertex(const std::vector<HRefitCand> &cands);
    std::vector<HRefitCand> UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &VertexPos);

    void setLearningRate(double val) { fLearningRate = val; }
    void setNumberOfIterations(int val) { fNumIterations = val; }

    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }
    TVector3 getVertex() const { return fVertex; } // Function that the user should use in the analysis macro

    // The functions below are functions to obtain information about distances
    // from the analysis macro
    double getDistanceBetweenFittedParticles() const { return fDistanceParticleToParticle; }
    double getDistanceFirstParticleVertex() const { return fDistParticle1Vertex; }
    double getDistanceSecondParticleVertex() const { return fDistParticle2Vertex; }
    double getDistanceFirstParticleOrigin() const { return fDistParticle1Origin; }
    double getDistanceSecondParticleOrigin() const { return fDistParticle2Origin; }

    void setLambdaCandidate(double valTheta, double valPhi, double valR, double ValZ);
    HVirtualCand getLambdaCandidate() { return fLambdaCandidate; }

    bool isConverged() const { return fConverged; }
    int getIteration() const { return fIteration; }
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }
    bool fit();
    void setVerbosity(int val) { fVerbose = val; }

    void setPhiOriginal(double val1, double val2){
        fPhi1Original=val1;
        fPhi2Original=val2;
        }

    double fDistParticle1Vertex;
    double fDistParticle2Vertex;
    double fDistParticle1Origin;
    double fDistParticle2Origin;

    // J.R. The following line was present in the
    // original code. It does not seen to work with the
    // HYDRA version needed to run this code
    //[[deprecated]]

    HRefitCand getDaughter(int val);

    void update();

protected:
    void updateDaughters();
};

#endif /* HVERTEXFITTER_H */
