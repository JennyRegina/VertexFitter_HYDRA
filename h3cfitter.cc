#include "h3cfitter.h"

const size_t cov_dim = 5;

H3cFitter::H3cFitter(const std::vector<HRefitCand>& cands, HRefitCand& mother) : 
    fCands(cands),
    fMother(mother)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = (fN+1) * cov_dim -1;       //Dimension of covariance matrix (number of measured variables). Mother momentum is not measured

    y.ResizeTo(fyDim, 1); 
    p.ResizeTo(1,1);
    V.ResizeTo(fyDim, fyDim);
    Vp.ResizeTo(1,1);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;

    // set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) //for daughters
    {
        HRefitCand cand = cands[ix];

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
	
	//for mother
    HRefitCand cand = fMother;

    y(fN * cov_dim, 0) = cand.Theta();
    y(1 + fN * cov_dim, 0) = cand.Phi();
    y(2 + fN * cov_dim, 0) = cand.getR();
    y(3 + fN * cov_dim, 0) = cand.getZ();
    fM.push_back(mother.M());

    TMatrixD covariance = mother.getCovariance();
    V(0 + fN * cov_dim, 0 + fN * cov_dim) = covariance(1, 1);
    V(1 + fN * cov_dim, 1 + fN * cov_dim) = covariance(2, 2);
    V(2 + fN * cov_dim, 2 + fN * cov_dim) = covariance(3, 3);
    V(3 + fN * cov_dim, 3 + fN * cov_dim) = covariance(4, 4);
}

void H3cFitter::add3Constraint()
{
    fNdf += 3;
    f3Constraint = true;
}

//Calculate starting value of mother 1/p from energy conservation constraint
TMatrixD H3cFitter::calcMotherMom(const TMatrixD& m_iter)
{
    TMatrix xi(1,1);

    for(int q=0; q<fN; q++){
        xi(0, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
    }
    xi(0, 0) = pow(xi(0,0),2)-pow(fM[fN],2);

    xi(0, 0) = 1./sqrt(xi(0,0));

    return xi;
}

//Constraint equations
TMatrixD H3cFitter::f_eval(const TMatrixD& m_iter, const TMatrixD& xi_iter)
{
    TMatrixD d;

    if (f3Constraint)
    {
	//mother
	d.ResizeTo(4, 1);
        d(0,0) = -1. / xi_iter(0, 0)*sin(m_iter(0 + fN, 0))*cos(m_iter(1 + fN, 0));
        d(1,0) = -1. / xi_iter(0, 0)*sin(m_iter(0 + fN, 0))*sin(m_iter(1 + fN, 0));
        d(2,0) = -1. / xi_iter(0, 0)*cos(m_iter(0 + fN, 0));
        d(3,0) = -sqrt(pow((1. / xi_iter(0, 0)), 2) + pow(fM[fN], 2));
        
	//daughters
        for(int q=0; q<fN; q++){
            d(0,0) += 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*cos(m_iter(2 + q * cov_dim, 0));
            d(1,0) += 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*sin(m_iter(2 + q * cov_dim, 0));
            d(2,0) += 1. / m_iter(0 + q * cov_dim, 0)*cos(m_iter(1 + q * cov_dim, 0));
            d(3,0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
    }

    return d;
}

//Jacobian (derivative of constraint equations with respect to measured variables)
TMatrixD H3cFitter::Feta_eval(const TMatrixD& m_iter, const TMatrixD& xi_iter)
{

    TMatrixD H;

    if (f3Constraint)
    {
        H.ResizeTo(4, fyDim);
        H.Zero();
		
		//Daughter variables
        for(int q=0; q<fN; q++){
			//d(1/p)
            H(0, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),2) * cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1./pow(m_iter(0 + q * cov_dim, 0),3) * 1./ sqrt(pow(1./(m_iter(0 + q * cov_dim, 0)),2) + pow(fM[q], 2));

            //dtht
            H(0, 1 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0));

            //dphi
            H(0, 2 + q * cov_dim) = -1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1./m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
        }
		
		//Mother variables
		//dtht
        H(0, fN * cov_dim) = -1./xi_iter(0, 0) * cos(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
        H(1, fN * cov_dim) = -1./xi_iter(0, 0) * cos(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(2, fN * cov_dim) = 1./xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0));
        //dphi
        H(0, 1 + fN * cov_dim) = 1./xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(1, 1 + fN * cov_dim) = -1./xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
    }

    return H;
}

//Jacobian (derivative of constraint equations with respect to unmeasured variables)
TMatrixD H3cFitter::Fxi_eval(const TMatrixD& m_iter, const TMatrixD& xi_iter)
{

    TMatrixD H;

    if (f3Constraint)
    {
        H.ResizeTo(4, 1);
        H.Zero();

        //d(1/p)
        H(0, 0) = 1./pow(xi_iter(0, 0),2) * sin(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
        H(1, 0) = 1./pow(xi_iter(0, 0),2) * sin(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(2, 0) = 1./pow(xi_iter(0, 0),2) * cos(m_iter(0 + fN * cov_dim, 0));
        H(3, 0) = 1./pow(xi_iter(0, 0),3) * 1./ sqrt(pow(m_iter(0, 0),2) + pow(fM[fN], 2));
    }

    return H;
}

bool H3cFitter::fit(double lr, Int_t maxItr)
{
   // double lr = 0.5;
    TMatrixD alpha0(fyDim, 1), alpha(fyDim, 1);
    TMatrixD xi0(1, 1), xi(1, 1);
    TMatrixD A0(y), V0(V);
    
    TMatrixD V0_inv(V);
    V0_inv.Invert();

    alpha0 = y;
    alpha = alpha0;

    //Calculate unknown mother momentum from constraints
    xi0 = calcMotherMom(alpha0);
    xi = xi0;
    double chi2 = 1e6;
    //cout << " calc f " << endl;
    TMatrixD d = f_eval(alpha, xi);
    //cout << " calc Feta" << endl;
    TMatrixD D = Feta_eval(alpha, xi);
    //cout << " calc Fxi" << endl;
    TMatrixD D_xi = Fxi_eval(alpha, xi);
    //cout << " start fitting " << endl;

    for (int q = 0; q < maxItr; q++)
    {   
        TMatrixD delta_alpha = alpha0 - alpha;
        //calc r
        TMatrixD r = d + D * delta_alpha;
        TMatrixD DT(D.GetNcols(), D.GetNrows());
        DT.Transpose(D);
        TMatrixD DT_xi(D_xi.GetNcols(), D_xi.GetNrows());
        DT_xi.Transpose(D_xi);
	    //calc S
        TMatrixD VD = D * V0 * DT;
        VD.Invert();
        TMatrixD VDD = DT_xi * VD * D_xi;
        VDD.Invert();

        //calculate values for next iteration
        TMatrixD neu_xi = xi - lr * VDD * DT_xi * VD * r;
        TMatrixD delta_xi = neu_xi - xi;
	//cout << " calc lambda " << endl;
        TMatrixD lambda = VD * (r + D_xi * delta_xi);
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fyDim, 1);
        neu_alpha = alpha0 - lr * V0 * DT * lambda; //oder alpha?

        //Update covariance - could be done after fitting
        TMatrixD matrix = DT*VD*D_xi;
        TMatrixD matrixT(matrix.GetNcols(), matrix.GetNrows());
        matrixT.Transpose(matrix);
        TMatrixD invertedMatrix = DT_xi*VD*D_xi;
        invertedMatrix.Invert();
        V = V0 - lr * V0 * (DT * VD * D - (matrix*invertedMatrix*matrixT)) * V0;
        Vp = invertedMatrix;

        //Calculate new chi2
        TMatrixD chisqrd(1,1);
        TMatrixD delta_alphaT(delta_alpha.GetNcols(), delta_alpha.GetNrows());
        delta_alphaT.Transpose(delta_alpha);
        TMatrixD two(1,1);
        two(0,0) = 2;
        chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * d;


        // for checking convergence
        // three parameters are checked
        // 1. difference between measurements (for successive iterations) y
        // 2. difference between constraints (for successive iterations)  d
        // 3. difference between chi2 (for successive iterations)  chisqrd
        // check converge for 'y' measurements
        /*
        double sum0 = 0;
        for(uint p=0; p<(fN*cov_dim); p++){
            sum0 += (neu_alpha(p,0)-alpha(p,0))*(neu_alpha(p,0)-alpha(p,0));
        }

        double d_const = fabs(d(0,0));
        */
        if(fabs(chi2-chisqrd(0,0))<1){
            fIteration = q;
            fConverged = true;
            chi2 = chisqrd(0,0);
            alpha = neu_alpha;
            xi = neu_xi;
            //V = V - lr * V * DT * VD * D * V;
	    //cout << q << endl;
            break;
        }
        
	
        chi2 = chisqrd(0,0);
        alpha = neu_alpha; 
        xi = neu_xi;
        D = Feta_eval(alpha,xi);
        D_xi = Fxi_eval(alpha, xi);
        d = f_eval(alpha,xi);
    }

    y = alpha;
    p = xi;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(fyDim, fyDim);
    for (int b = 0; b < (fyDim); b++)
        fPull(b, b) = -10000;

    if (true)
    {
        for (int b = 0; b < (fyDim); b++)
        {
            double num = A0(b, 0) - alpha(b, 0);
            double dem = V0(b, b) - V(b, b);
            if (dem > 0) { fPull(b, b) = num / std::sqrt(dem); }
        }
    }

    updateDaughters();
    updateMother();

    //return fConverged; // for number of iterations greater than 1
    return true; // for number of iterations equal to 1
}

HRefitCand H3cFitter::getDaughter(int val)
{
    return fCands[val];
}

HRefitCand H3cFitter::getMother()
{
    return fMother;
}

void H3cFitter::updateDaughters()
{
    for (int val = 0; val < fN; ++val)
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
        cand.SetXYZM(Px, Py, Pz, M);
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

void H3cFitter::updateMother()
{
    HRefitCand& mother = fMother;
    double Px = (1. / p(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::cos(y(1 + fN * cov_dim, 0));
    double Py = (1. / p(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::sin(y(1 + fN * cov_dim, 0));
    double Pz =
        (1. / p(0, 0)) * std::cos(y(0 + fN * cov_dim, 0));
    double M = fM[fN];
    mother.SetXYZM(Px, Py, Pz, M);
    mother.setR(y(2 + fN * cov_dim, 0));
    mother.setZ(y(3 + fN * cov_dim, 0));

    // ---------------------------------------------------------------------------
    // set covariance
    // ---------------------------------------------------------------------------
    TMatrixD cov(5, 5);
    cov(0, 0) = Vp(0, 0);
    cov(1, 1) = V(1 + fN * cov_dim, 1 + fN * cov_dim);
    cov(2, 2) = V(2 + fN * cov_dim, 2 + fN * cov_dim);
    cov(3, 3) = V(3 + fN * cov_dim, 3 + fN * cov_dim);
    cov(4, 4) = V(4 + fN * cov_dim, 4 + fN * cov_dim);
    mother.setCovariance(cov);
    // ---------------------------------------------------------------------------
}

void H3cFitter::update()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand& cand = fCands[val];
        cand.update();
    }
}
