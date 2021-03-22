#include "hvertexfitter.h"

const size_t cov_dim = 5;

HVertexFitter::HVertexFitter(const std::vector<HRefitCand> &cands) : 
    fCands(cands), 
    fVerbose(0), 
    fLearningRate(1), 
    fNumIterations(10)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = fN * cov_dim;       //Dimension of full covariance matrix (number of measured variables x cov_dim)

    y.ResizeTo(fyDim, 1);
    V.ResizeTo(fyDim, fyDim);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    fVtxConstraint = false;

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix = 0; ix < fN; ix++)
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
}

HVertexFitter::HVertexFitter(const std::vector<HRefitCand>& cands, HRefitCand& mother) : 
    fCands(cands),
    fMother(mother),
    fVerbose(0), 
    fLearningRate(1), 
    fNumIterations(10)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = (fN+1) * cov_dim -1;       //Dimension of full covariance matrix (number of measured variables x cov_dim). Mother momentum is not measured

    y.ResizeTo(fyDim, 1); 
    x.ResizeTo(1,1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(1,1);

    y.Zero();
    V.Zero();
    x.Zero();
    Vx.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    fVtxConstraint = false;

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
    fM.push_back(cand.M());

    TMatrixD covariance = cand.getCovariance();
    V(0 + fN * cov_dim, 0 + fN * cov_dim) = covariance(1, 1);
    V(1 + fN * cov_dim, 1 + fN * cov_dim) = covariance(2, 2);
    V(2 + fN * cov_dim, 2 + fN * cov_dim) = covariance(3, 3);
    V(3 + fN * cov_dim, 3 + fN * cov_dim) = covariance(4, 4);
}


void HVertexFitter::add3Constraint()
{
    fNdf += 3;
    f3Constraint = true;
}

void HVertexFitter::addVertexConstraint()
{
    fNdf += 1;
    fVtxConstraint = true;
}

//Constraint equations
TMatrixD HVertexFitter::f_eval(const TMatrixD &m_iter, const TMatrixD& xi_iter)
{
    TMatrixD d;

    if(fVtxConstraint){
        d.ResizeTo(1, 1);
        TVector3 base_1, base_2, dir_1, dir_2;

        // J.R The cobe block below is the original base vectors where they are calculated from the track parameters
        base_1.SetXYZ(
            m_iter(3 + 0 * cov_dim, 0) *
                std::cos(m_iter(2 + 0 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(3 + 0 * cov_dim, 0) *
                std::sin(m_iter(2 + 0 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(4 + 0 * cov_dim, 0));
        base_2.SetXYZ(
            m_iter(3 + 1 * cov_dim, 0) *
                std::cos(m_iter(2 + 1 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(3 + 1 * cov_dim, 0) *
                std::sin(m_iter(2 + 1 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(4 + 1 * cov_dim, 0));

        // The base vectors below are calculated from two constant phi,
        // These are the ones for the original base vector pointing from the origin
        // to the point along the track closest to the beam line

        /*     base_1.SetXYZ(
            m_iter(3 + 0 * cov_dim, 0) *
                std::cos(fPhi1Original + TMath::PiOver2()),
            m_iter(3 + 0 * cov_dim, 0) *
                std::sin(fPhi1Original + TMath::PiOver2()),
            m_iter(4 + 0 * cov_dim, 0));
        base_2.SetXYZ(
            m_iter(3 + 1 * cov_dim, 0) *
                std::cos(fPhi2Original + TMath::PiOver2()),
            m_iter(3 + 1 * cov_dim, 0) *
                std::sin(fPhi2Original + TMath::PiOver2()),
            m_iter(4 + 1 * cov_dim, 0)); */

        dir_1.SetXYZ(std::sin(m_iter(1 + 0 * cov_dim, 0)) *
                        std::cos(m_iter(2 + 0 * cov_dim, 0)),
                    std::sin(m_iter(1 + 0 * cov_dim, 0)) *
                        std::sin(m_iter(2 + 0 * cov_dim, 0)),
                    std::cos(m_iter(1 + 0 * cov_dim, 0)));
        dir_2.SetXYZ(std::sin(m_iter(1 + 1 * cov_dim, 0)) *
                        std::cos(m_iter(2 + 1 * cov_dim, 0)),
                    std::sin(m_iter(1 + 1 * cov_dim, 0)) *
                        std::sin(m_iter(2 + 1 * cov_dim, 0)),
                    std::cos(m_iter(1 + 1 * cov_dim, 0)));

        d(0, 0) = std::fabs((dir_1.Cross(dir_2)).Dot((base_1 - base_2)));
    }
    
    if (f3Constraint)
    {
	//mother
	d.ResizeTo(4, 1);
        d(0,0) = -1. / xi_iter(0, 0)*sin(m_iter(0 + fN * cov_dim, 0))*cos(m_iter(1 + fN * cov_dim, 0));
        d(1,0) = -1. / xi_iter(0, 0)*sin(m_iter(0 + fN * cov_dim, 0))*sin(m_iter(1 + fN * cov_dim, 0));
        d(2,0) = -1. / xi_iter(0, 0)*cos(m_iter(0 + fN * cov_dim, 0));
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
TMatrixD HVertexFitter::Feta_eval(const TMatrixD &m_iter, const TMatrixD& xi_iter)
{
    TMatrixD H;

    if(fVtxConstraint){
        H.ResizeTo(1, fyDim);
        H.Zero();

        H(0, 0) = 0.;
        H(0, 5) = 0.;

        H(0, 1) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::sin(m_iter(2, 0)) -
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) -
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::sin(m_iter(1, 0)) -
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) -
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::sin(m_iter(1, 0)) +
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::sin(m_iter(1, 0)) +
                (m_iter(4, 0) - m_iter(9, 0)) *
                    (std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                    std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 6) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) -
                m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) + //
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) +
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) +
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) +
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) -
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) -
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) +
                (m_iter(4, 0) - m_iter(9, 0)) *
                    (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 2) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) -
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) -
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) -
                m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) - //
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                (m_iter(9, 0) - m_iter(4, 0)) *
                    (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 7) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) +
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) -
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) -
                m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) +
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::cos(m_iter(6, 0)) +
                m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) -
                m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                    std::cos(m_iter(1, 0)) +
                (m_iter(4, 0) - m_iter(9, 0)) *
                    (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) -
                    std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)));

        H(0, 3) = std::cos(m_iter(2, 0) + pi2) *
                    (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) -
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                        std::cos(m_iter(1, 0))) -
                std::sin(m_iter(2, 0) + pi2) *
                    (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) -
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                        std::cos(m_iter(1, 0)));

        H(0, 8) = -1 * std::cos(m_iter(7, 0) + pi2) *
                    (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) -
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                        std::cos(m_iter(1, 0))) +
                std::sin(m_iter(7, 0) + pi2) *
                    (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                        std::cos(m_iter(6, 0)) -
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                        std::cos(m_iter(1, 0)));

        H(0, 4) = std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));

        H(0, 9) = -1 * std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) +
                std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                    std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));

    }

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
TMatrixD HVertexFitter::Fxi_eval(const TMatrixD& m_iter, const TMatrixD& xi_iter)
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

bool HVertexFitter::fit()
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFitter::fit() -----------" << std::endl;
        std::cout << "Vertex constraint set: " << fVtxConstraint << std::endl;
        std::cout << "3C set: " << f3Constraint << std::endl;
        std::cout << "" << std::endl;
    }

    double lr = fLearningRate;
    TMatrixD alpha0(fyDim, 1), alpha(fyDim, 1);
    TMatrixD xi0(1, 1), xi(1, 1), neu_xi(1, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;

    // Calculating the inverse of the original covariance matrix that is not changed in the iterations
    TMatrixD V0_inv(V); // J.R New
    V0_inv.Invert();    // J.R New

    xi0.Zero();
    xi.Zero();
    neu_xi.Zero();

    if(f3Constraint){
        xi0(0,0) = 1/fMother.P();
        xi = xi0;
    }

    double chi2 = 1e6;
    TMatrixD D = Feta_eval(alpha, xi);
    TMatrixD DT(D.GetNcols(), D.GetNrows());
    TMatrixD d = f_eval(alpha, xi);
    TMatrixD D_xi(d.GetNrows(), 1), DT_xi(1, d.GetNrows());       //check dimension if other fitters are added
    D_xi.Zero();
    DT_xi.Zero();
    if(f3Constraint) D_xi = Fxi_eval(alpha, xi);
    TMatrixD VD(D.GetNrows(), D.GetNrows());
    VD.Zero();
    TMatrixD VDD(D_xi.GetNrows(), D_xi.GetNrows());
    VDD.Zero();

    for (int q = 0; q < fNumIterations; q++)
    {
        TMatrixD delta_alpha = alpha0 - alpha;
        //calc r
        TMatrixD r = d + D * delta_alpha;
        DT.Transpose(D);
        //calc S
        VD = D * V0 * DT;
        VD.Invert();
        if(f3Constraint){
            DT_xi.Transpose(D_xi);
            VDD = DT_xi * VD * D_xi;
            VDD.Invert();
        }

        //calculate values for next iteration
        TMatrixD lambda(d);       //Lagrange multiplier
        lambda.Zero();
        if(f3Constraint){
            neu_xi = xi - lr * VDD * DT_xi * VD * r;
        }
        TMatrixD delta_xi = neu_xi - xi;
        lambda = VD * (r + D_xi * delta_xi);
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fyDim, 1);
        neu_alpha = alpha0 - lr * V0 * DT * lambda; 
        
        //Calculate new chi2
        TMatrixD chisqrd(1, 1);
        TMatrixD delta_alphaT(delta_alpha.GetNcols(), delta_alpha.GetNrows());
        delta_alphaT.Transpose(delta_alpha);
        TMatrixD two(1, 1);
        two(0, 0) = 2;
        chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * d;
        
        //for checking convergence
        if(fabs(chi2-chisqrd(0,0))<1){
            fIteration = q;
            fConverged = true;
            chi2 = chisqrd(0,0);
            alpha = neu_alpha;
            if(f3Constraint) xi = neu_xi;
            break;
        }

        if (fVerbose > 0)
        {
            std::cout << "Iteration: " << q << std::endl;
            std::cout << "Printing d: " << std::endl; 
            d.Print();
        }

        fIteration = q;
        chi2 = chisqrd(0, 0);
        alpha = neu_alpha;
        if(f3Constraint) xi = neu_xi;
        D = Feta_eval(alpha, xi);
        if(f3Constraint) D_xi = Fxi_eval(alpha, xi);
        d = f_eval(alpha, xi);
    }

    y = alpha;
    x = xi;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    //Update covariance
    if(f3Constraint){
        TMatrixD matrix = DT*VD*D_xi;
        TMatrixD matrixT(matrix.GetNcols(), matrix.GetNrows());
        matrixT.Transpose(matrix);
        TMatrixD invertedMatrix = DT_xi*VD*D_xi;
        invertedMatrix.Invert();
        V = V0 - lr * V0 * (DT * VD * D - (matrix*invertedMatrix*matrixT)) * V0;
        Vx = invertedMatrix;
    }
    if(fVtxConstraint) V = V0 - lr * V0 * DT * VD * D * V0;

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
            if (dem > 0)
            {
                fPull(b, b) = num / std::sqrt(dem);
            }
        }
    }

    updateDaughters();
    if(f3Constraint) updateMother();

    // return fConverged; // for number of iterations greater than 1
    return true; // for number of iterations equal to 1
}

HRefitCand HVertexFitter::getDaughter(int val)
{
    return fCands[val];
}

HRefitCand HVertexFitter::getMother()
{
    return fMother;
}

void HVertexFitter::updateDaughters()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand &cand = fCands[val];
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

void HVertexFitter::updateMother()
{
    HRefitCand& mother = fMother;
    double Px = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::cos(y(1 + fN * cov_dim, 0));
    double Py = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::sin(y(1 + fN * cov_dim, 0));
    double Pz =
        (1. / x(0, 0)) * std::cos(y(0 + fN * cov_dim, 0));
    double M = fM[fN];
    mother.SetXYZM(Px, Py, Pz, M);
    mother.setR(y(2 + fN * cov_dim, 0));
    mother.setZ(y(3 + fN * cov_dim, 0));

    // ---------------------------------------------------------------------------
    // set covariance
    // ---------------------------------------------------------------------------
    TMatrixD cov(5, 5);
    cov(0, 0) = Vx(0, 0);
    cov(1, 1) = V(1 + fN * cov_dim, 1 + fN * cov_dim);
    cov(2, 2) = V(2 + fN * cov_dim, 2 + fN * cov_dim);
    cov(3, 3) = V(3 + fN * cov_dim, 3 + fN * cov_dim);
    cov(4, 4) = V(4 + fN * cov_dim, 4 + fN * cov_dim);
    mother.setCovariance(cov);
    // ---------------------------------------------------------------------------
}

void HVertexFitter::update()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand &cand = fCands[val];
        cand.update();
    }
}
