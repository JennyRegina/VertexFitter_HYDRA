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
#include "hgeomvector.h" // J.R added for vertex calculation
#include "hparticletool.h" // J.R added for vertex calculation

using std::cout;
using std::endl;

const double pi2 = TMath::PiOver2();

template<typename T>
void Print(T const &matrix)
{
    Int_t nrows = matrix.GetNrows();
    Int_t ncols = matrix.GetNcols();
    
    cout << "shape(" << nrows << "," << ncols << ")" << endl;

    for (Int_t i=0; i<nrows; i++){
        for (Int_t j=0; j<ncols; j++){
            Double_t element = matrix(i,j);
            if ( TMath::Abs(element) < 1e-10 ) element=0.;
            if ( element >= 0.)
              cout << " " << std::fixed << std::setw(8) << std::scientific << element << " ";
            else
              cout << std::fixed << std::setw(8) << std::scientific << element << " ";
        }
        cout << endl;
    }

    cout << endl;
}

class HVertexFitter {
private:
    TMatrixD y, V, fPull;
    double fChi2, fProb;
    bool   fConverged;
    int    fIteration, fN;
    std::vector<HRefitCand> fCands;

    // data members for constraints
    double fMass;
    int fNdf;
    std::vector<double> fM;
    TLorentzVector fInit;
    bool fVtxConstraint;  
    double fVtxPos;
    TVector3 fVertex;
    

public:
    HVertexFitter(const std::vector<HRefitCand> & cands);
    ~HVertexFitter(){};
    TMatrixD f_eval(const TMatrixD &m_iter);
    TMatrixD Feta_eval(const TMatrixD &miter);
    // The addVtxConstraint function should work in two different ways:
    // 1. Take the fitted vertex and adjust track parameters to be coming from there
    // 2. Possibility: adjust the reconstructed vertex position  
    void   addVtxConstraint();
    TVector3 findVertex();
    //TVector3 findVertex(const std::vector<HRefitCand> & cands) {return fVertex;} // Function to calculate and return vertex
    double getChi2() const {return fChi2;}
    double getProb() const {return fProb;}
    double getPull(int val=0){return fPull(val,val);}
    bool   isConverged() const {return fConverged;}
    int    getIteration() const {return fIteration;}
    void   setCovariance(TMatrixD &val){V=val;}
    void   setMeasurement(TMatrixD &val){y=val;}

    bool fit();

    //[[deprecated]]
    HRefitCand getDaughter(int val);

    void update();

protected:
    void updateDaughters();
};

#endif /* HVERTEXFITTER_H */
