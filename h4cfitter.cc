#include "h4cfitter.h"

const size_t cov_dim = 5;

H4cFitter::H4cFitter(const std::vector<HRefitCand>& cands, HRefitCand& mother) : 
    fCands(cands)
   // fMother(mother)
{
    // fNdau is the number of daughters e.g. (L->ppi-) n=2
    fNdau = cands.size();
    fNdau++;
    fCands.push_back(mother);
    fWiggleMoth = true;

    y.ResizeTo(fNdau * cov_dim, 1);
    V.ResizeTo(fNdau * cov_dim, fNdau * cov_dim);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    fLv4C = TLorentzVector(0, 0, 0, 0);

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix = 0; ix < fNdau; ix++)
    {
        HRefitCand cand = cands[ix];
        //fM[ix]=cand.M();

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        if(ix==fNdau-1){
            y(1 + ix * cov_dim, 0) = pi + cand.Theta();
        }else{
            y(1 + ix * cov_dim, 0) = cand.Theta();
        }
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }

    f4Constraint = false;
}

//H4cFitter::H4cFitter(const std::vector<HRefitCand>& cands) : 
//    fCands(cands)
H4cFitter::H4cFitter(const std::vector<HRefitCand>& cands, TLorentzVector& lv) : 
    fCands(cands)
{
    //fLv4C = TLorentzVector(0,0,4337.96,2*938.272+3500) ;
    fLv4C = lv;
    
    // fNdau is the number of daughters e.g. (L->ppi-) n=2
    fNdau = cands.size();
    fWiggleMoth = false;

    y.ResizeTo(fNdau * cov_dim, 1);
    V.ResizeTo(fNdau * cov_dim, fNdau * cov_dim);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix = 0; ix < fNdau; ix++)
    {
        HRefitCand cand = cands[ix];
        //fM[ix]=cand.M();
        
        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }

    f4Constraint = false;
}

void H4cFitter::add4Constraint()
{
    fNdf += 4;
    f4Constraint = true;
}

TMatrixD H4cFitter::f_eval(const TMatrixD& m_iter)
{
    TMatrixD d;

    // 4 constraint
    if (f4Constraint)
    {
	
	d.ResizeTo(fNdf, 1); //(4,1)
        d(0, 0) = -fLv4C.Px();
        d(1, 0) = -fLv4C.Py();
        d(2, 0) = -fLv4C.Pz();
        d(3, 0) = -fLv4C.E();
	
        for(int q=0; q<fNdau; q++){
            d(0,0) += 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*cos(m_iter(2 + q * cov_dim, 0));
            d(1,0) += 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*sin(m_iter(2 + q * cov_dim, 0));
            d(2,0) += 1. / m_iter(0 + q * cov_dim, 0)*cos(m_iter(1 + q * cov_dim, 0));
            d(3,0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }

    }

    return d;
}

TMatrixD H4cFitter::Feta_eval(const TMatrixD& m_iter)
{

    TMatrixD H;

    // 4constraint
    if (f4Constraint)
    {
        H.ResizeTo(4, fNdau * cov_dim);
        H.Zero();

        for(int q=0; q<fNdau; q++){
            H(0, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),3) * 1./ sqrt(pow(1./(m_iter(0 + q * cov_dim, 0)),2) + pow(fM[q], 2));

            H(0, 1 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0));

            H(0, 2 + q * cov_dim) = -1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
        }
    }

    return H;
}

bool H4cFitter::fit(double lr, Int_t maxItr)
{
   // double lr = 0.5;
    TMatrixD alpha0(fNdau * cov_dim, 1), alpha(fNdau * cov_dim, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    double chi2 = 1e6;
    //cout << " calc Feta" << endl;
    TMatrixD D = Feta_eval(alpha);
    //cout << " calc f " << endl;
    TMatrixD d = f_eval(alpha);
    //cout << " start fitting " << endl;

    for (int q = 0; q < maxItr; q++)
    {
        TMatrixD DT(D.GetNcols(), D.GetNrows());
        DT.Transpose(D);
	//cout << " calc D " << endl;
        TMatrixD VD = D * V * DT;
        VD.Invert();

        //TMatrixD delta_alpha = alpha - alpha0;
        TMatrixD delta_alpha = y - alpha;
	//cout << " calc lambda " << endl;
        TMatrixD lambda = VD * D * delta_alpha + VD * d; 
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fNdau * cov_dim, 1);
        neu_alpha = alpha - lr * V * DT * lambda;

        double chisqrd = 0.;

        for (int p = 0; p < lambda.GetNrows(); p++)
        {
            chisqrd += lambdaT(0, p) * d(p, 0);
        }

        // for checking convergence
        // three parameters are checked
        // 1. difference between measurements (for successive iterations) y
        // 2. difference between constraints (for successive iterations)  d
        // 3. difference between chi2 (for successive iterations)  chisqrd
        // check converge for 'y' measurements
        /*
        double sum0 = 0;
        for(uint p=0; p<(fNdau*cov_dim); p++){
            sum0 += (neu_alpha(p,0)-alpha(p,0))*(neu_alpha(p,0)-alpha(p,0));
        }

        double d_const = fabs(d(0,0));
        */
        if(fabs(chi2-chisqrd)<1){
            fIteration = q;
            fConverged = true;
            chi2 = chisqrd;
            alpha0 = alpha;
            alpha = neu_alpha;
            V = V - lr * V * DT * VD * D * V;
	    //cout << q << endl;
            break;
        }
        
	
        chi2 = chisqrd;
        alpha0 = alpha;
        alpha = neu_alpha;
        V = V - lr * V * DT * VD * D * V;
        D = Feta_eval(alpha);
        d = f_eval(alpha);
    }

    y = alpha;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(fNdau * cov_dim, fNdau * cov_dim);
    for (uint b = 0; b < (fNdau * cov_dim); b++)
        fPull(b, b) = -10000;

    if (true)
    {
        for (uint b = 0; b < (fNdau * cov_dim); b++)
        {
            double num = A0(b, 0) - alpha(b, 0);
            double dem = V0(b, b) - V(b, b);
            if (dem > 0) { fPull(b, b) = num / std::sqrt(dem); }
        }
    }

    updateDaughters();

    //return fConverged; // for number of iterations greater than 1
    return true; // for number of iterations equal to 1
}

HRefitCand H4cFitter::getDaughter(int val)
{
    return fCands[val];
}

void H4cFitter::updateDaughters()
{
    for (int val = 0; val < fNdau; ++val)
    {
        HRefitCand& cand = fCands[val];
        double Px = (1. / y(0 + val * cov_dim, 0)) *
                    std::sin(y(1 + val * cov_dim, 0)) *
                    std::cos(y(2 + val * cov_dim, 0));
        double Py = (1. / y(0 + val * cov_dim, 0)) *
                    std::sin(y(1 + val * cov_dim, 0)) *
                    std::sin(y(2 + val * cov_dim, 0));
        double Pz =
            (1. / y(0 + val * cov_dim, 0)) * std::cos(y(1 + val * cov_dim, 0));
        double M = fM[val];
        if(fWiggleMoth && val==fNdau-1){
            cand.SetXYZM(-Px, -Py, -Pz, M);
        }else{
            cand.SetXYZM(Px, Py, Pz, M);
        }
        cand.setR(y(3 + val * cov_dim, 0));
        cand.setZ(y(4 + val * cov_dim, 0));

        // ---------------------------------------------------------------------------
        // set covariance
        // ---------------------------------------------------------------------------
        TMatrixD cov(5, 5);
        cov(0, 0) = V(0 + val * cov_dim, 0 + val * cov_dim);
        cov(1, 1) = V(1 + val * cov_dim, 1 + val * cov_dim);
        cov(2, 2) = V(2 + val * cov_dim, 2 + val * cov_dim);
        cov(3, 3) = V(3 + val * cov_dim, 3 + val * cov_dim);
        cov(4, 4) = V(4 + val * cov_dim, 4 + val * cov_dim);
        cand.setCovariance(cov);
        // ---------------------------------------------------------------------------
    }
}

void H4cFitter::update()
{
    for (int val = 0; val < fNdau; ++val)
    {
        HRefitCand& cand = fCands[val];
        cand.update();
    }
}
