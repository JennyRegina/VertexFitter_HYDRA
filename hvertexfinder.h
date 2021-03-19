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

using std::cout;
using std::endl;

//const double pi2 = TMath::PiOver2();
/*
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
*/
class HVertexFinder
{
private:

    std::vector<HRefitCand> fCands;

    // data members for constraints
    double fMass;
    int fNdf;

    std::vector<double> fM;
    TLorentzVector fInit;
    bool fVtxConstraint;
    double fVtxPos;
    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    int fVerbose;

    double fDistanceParticleToParticle;
    double fDistanceParticleToVertex;

    double fPhi1Original;
    double fPhi2Original;
    HVirtualCand fLambdaCandidate;

public:
    HVertexFinder(const std::vector<HRefitCand> &cands);
    ~HVertexFinder(){};
    TVector3 findVertex(const std::vector<HRefitCand> &cands);
    std::vector<HRefitCand> UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &VertexPos);

    TVector3 getVertex() const { return fVertex; } // Function that the user should use in the analysis macro

    // The functions below are functions to obtain information about distances
    // from the analysis macro
    double getDistanceBetweenFittedParticles() const { return fDistanceParticleToParticle; }
    double getDistanceFirstParticleVertex() const { return fDistParticle1Vertex; }
    double getDistanceSecondParticleVertex() const { return fDistParticle2Vertex; }
    double getDistanceFirstParticleOrigin() const { return fDistParticle1Origin; }
    double getDistanceSecondParticleOrigin() const { return fDistParticle2Origin; }

    // The first function is for creating a Lambda candidate if only information of the decay vertex is available
    // The second function is for creating the Lambda candidate if information about the primary vertex is also available

    void setLambdaCandidate(double valMomentum, double valTheta, double valPhi, double valR, double ValZ, TVector3 decayVertex);
    void  setLambdaCandidateFromPrimaryVtxInfo(double valMomentum, double valTheta, double valPhi, double valR, double ValZ, TVector3 primVtx);
    HVirtualCand getLambdaCandidate() { return fLambdaCandidate; }
    TMatrixD getCovarianceMatrixLambda() { return fCovarianceLambda; }

    void setVerbosity(int val) { fVerbose = val; }

    void setPhiOriginal(double val1, double val2){
        fPhi1Original=val1;
        fPhi2Original=val2;
        }

    double fDistParticle1Vertex;
    double fDistParticle2Vertex;
    double fDistParticle1Origin;
    double fDistParticle2Origin;

    //TMatrixD fCovarianceLambda(5, 5);
    TMatrixD fCovarianceLambda;

    // J.R. The following line was present in the
    // original code. It does not seen to work with the
    // HYDRA version needed to run this code
    //[[deprecated]]

    HRefitCand getDaughter(int val);

    // J.R. Variable noting if the primary vertex was found
    // This influences the calculation of the Lambda Candidate
    bool fPrimaryVertexFound;

protected:
    void updateDaughters();
};

#endif /* HVERTEXFINDER_H */
