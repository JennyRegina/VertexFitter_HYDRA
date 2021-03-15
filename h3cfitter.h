/**
 * H3CFitter.h
 *
 *
 */

#ifndef H3CFITTER_H
#define H3CFITTER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "TLorentzVector.h"

using std::cout;
using std::endl;

const double pi = TMath::Pi();

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

class H3cFitter {
private:
    TMatrixD y, V, fPull;
    double fChi2, fProb;
    bool   fConverged;
    int    fIteration, fNdau;
    std::vector<HRefitCand> fCands;

    // data members for constraints
    int fNdf;
    std::vector<double> fM;
    TLorentzVector fInit;
    TLorentzVector fLv4C;
    HRefitCand fMother;
    bool fWiggleMoth;
    bool f3Constraint;  
    //double fVtxPos;
    //TVector3 fVertex;
    

public:
    H3cFitter(const std::vector<HRefitCand> & cands, HRefitCand & mother);
    //H3cFitter(const std::vector<HRefitCand> & cands);
    ~H3cFitter(){};
    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD calcMotherMom(const TMatrixD& m_iter);
    void   add3Constraint();
    double getChi2() const {return fChi2;}
    double getProb() const {return fProb;}
    double getPull(int val=0){return fPull(val,val);}
    bool   isConverged() const {return fConverged;}
    int    getIteration() const {return fIteration;}
    void   setCovariance(TMatrixD &val){V=val;}
    void   setMeasurement(TMatrixD &val){y=val;}
    //void   appendCand(HRefitCand cand){fCands.push_back(cand);}

    bool fit(double lr, Int_t maxItr);

    HRefitCand getDaughter(int val);

    void update();

protected:
    void updateDaughters();
    void updateMother();
};

#endif /* H3CFITTER_H */
