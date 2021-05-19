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
    TMatrixD y, x, V, Vx, fPull;
    double fChi2, fProb;
    bool fConverged;
    int fIteration, fN, fyDim;    
    double fConvergenceCriteria;
    std::vector<HRefitCand> fCands;
    HRefitCand fMother;

    // data members for constraints
    int fNdf;
    std::vector<double> fM;
    TLorentzVector fLv4C;

    bool fVtxConstraint, f3Constraint, f4Constraint;
    int fVerbose;

    double fLearningRate;
    int fNumIterations;

public:
    HVertexFitter(const std::vector<HRefitCand> &cands);
    HVertexFitter(const std::vector<HRefitCand> &cands, HRefitCand &mother);
    ~HVertexFitter(){};

    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

    void add3Constraint();
    void add4Constraint();
    void addVertexConstraint();

    void setLearningRate(double val) { fLearningRate = val; }
    void setNumberOfIterations(int val) { fNumIterations = val; }
    void setConvergenceCriteria(double val) { fConvergenceCriteria = val; }

    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }

    bool isConverged() const { return fConverged; }
    int getIteration() const { return fIteration; }
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }

    bool fit();

    void setVerbosity(int val) { fVerbose = val; }

    HRefitCand getDaughter(int val);
    HRefitCand getMother();

    void update();

protected:
    void updateDaughters();
    void updateMother();
};

#endif /* HVERTEXFITTER_H */
