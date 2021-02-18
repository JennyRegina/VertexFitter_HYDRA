#include "hvertexfitter.h"

const size_t cov_dim = 5;

HVertexFitter::HVertexFitter(const std::vector<HRefitCand>& cands) : fCands(cands)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();

    y.ResizeTo(fN * cov_dim, 1);
    V.ResizeTo(fN * cov_dim, fN * cov_dim);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;

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

    fVtxConstraint = false;
}

// One will need to pass the correct objects to the function once it is decided what to use
TVector3 HVertexFitter::findVertex(const std::vector<HRefitCand> & cands){

//TVector3 HVertexFitter::findVertex(){

    // try to find the decay vertex from the most basic information so that 
    // the code is as independent as possible from objects used in the analysis

    double param_p_inv1, param_theta1, param_phi1, param_R1, param_Z1;
    double param_p_inv2, param_theta2, param_phi2, param_R2, param_Z2;

/*     for (int ix = 0; ix < fN; ix++)
    {
        HRefitCand cand = cands[ix];

        param_p_inv = 1. / cand.P();
        param_theta = cand.Theta();
        param_phi = cand.Phi();
        param_R = cand.getR();
        param_Z = cand.getZ();
    } */

        HRefitCand cand1 = cands[0];

        param_p_inv1 = 1. / cand.P();
        param_theta1 = cand.Theta();
        param_phi1 = cand.Phi();
        param_R1 = cand.getR();
        param_Z1 = cand.getZ();

        HRefitCand cand2 = cands[0];

        param_p_inv2 = 1. / cand.P();
        param_theta2 = cand.Theta();
        param_phi2 = cand.Phi();
        param_R2 = cand.getR();
        param_Z2 = cand.getZ();
    
    // Calculate the base and direction vectors of the two candidates 
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2(),
    param_R1 * std::sin(param_phi1 + TMath::PiOver2(),
    param_Z1);
    
    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2(),
    param_R2 * std::sin(param_phi2 + TMath::PiOver2(),
    param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1)*std::cos(param_theta1),
    std::sin(param_theta1)*std::sin(param_phi1),
    std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2)*std::cos(param_theta2),
    std::sin(param_theta2)*std::sin(param_phi2),
    std::cos(param_theta2));

    // Calculate the distance between the two tracks
    doouble dist = std::fabs((vtx_dir_1.Cross(vtx_dir_2)).Dot((vtx_base_1 - vtx_base_2)));
    // Keep the possibility to use this distance as a rough cut

    // Converting to HGeomVector in order to make use of built in funtions
    // NOTE: Base vectors will need to change for secondary vertex
    HGeomVector vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;
    
    vtx_geom_dir_1.setX(std::sin(param_theta1)*std::cos(param_theta1));
    vtx_geom_dir_1.setY(std::sin(param_theta1)*std::sin(param_phi1));
    vtx_geom_dir_1.setY(std::cos(param_theta1));
    vtx_geom_dir_2.setX(std::sin(param_theta2)*std::cos(param_theta2));
    vtx_geom_dir_2.setY(std::sin(param_theta2)*std::sin(param_phi2));
    vtx_geom_dir_2.setY(std::cos(param_theta2));

    vtx_geom_base_1.setX(param_R1 * std::cos(param_phi1 + TMath::PiOver2());
    vtx_geom_base_1.setX(param_R1 * std::sin(param_phi1 + TMath::PiOver2());
    vtx_geom_base_1.setX(param_Z1);

    vtx_geom_base_2.setX(param_R2 * std::cos(param_phi2 + TMath::PiOver2());
    vtx_geom_base_2.setX(param_R2 * std::sin(param_phi2 + TMath::PiOver2());
    vtx_geom_base_2.setX(param_Z2);

    HGeomVector vertex = HParticleTool::calculatePointOfClosestApproach(vtx_base_1, vtx_dir_1, vtx_dir_2, vtx_dir_2);

    fVertex.SetXYZ(vertex.X(),vertex.Y(),vertex.Z());

    return fVertex;
}

void HVertexFitter::addVtxConstraint()
{
    fNdf += 1;
    fVtxConstraint = true;
}

TMatrixD HVertexFitter::f_eval(const TMatrixD& m_iter)
{
    TMatrixD d;

    // vertex constraint
    if (fVtxConstraint)
    {
        d.ResizeTo(1, 1);
        TVector3 base_1, base_2, dir_1, dir_2;
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

    return d;
}

TMatrixD HVertexFitter::Feta_eval(const TMatrixD& m_iter)
{

    TMatrixD H;

    // vertex constraint
    if (fVtxConstraint)
    {
        H.ResizeTo(1, fN * cov_dim);
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

    return H;
}

bool HVertexFitter::fit()
{
    double lr = 1;
    TMatrixD alpha0(fN * cov_dim, 1), alpha(fN * cov_dim, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    double chi2 = 1e6;
    TMatrixD D = Feta_eval(alpha);
    TMatrixD d = f_eval(alpha);

    for (int q = 0; q < 5; q++)
    {
        TMatrixD DT(D.GetNcols(), D.GetNrows());
        DT.Transpose(D);
        TMatrixD VD = D * V * DT;
        VD.Invert();

        TMatrixD delta_alpha = alpha - alpha0;
        TMatrixD lambda = VD * D * delta_alpha + VD * d;
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fN * cov_dim, 1);
        neu_alpha = alpha - lr * V * DT * lambda;

        double chisqrd = 0.;

        for (int p = 0; p < lambda.GetNrows(); p++)
        {
            chisqrd = lambdaT(0, p) * d(p, 0);
        }

        /* for checking convergence
        // three parameters are checked
        // 1. difference between measurements (for successive iterations) y
        // 2. difference between constraints (for successive iterations)  d
        // 3. difference between chi2 (for successive iterations)  chisqrd
        // check converge for 'y' measurements
        double sum0 = 0;
        for(int p=0; p<(fN*5); p++){
            sum0 += (neu_alpha(p,0)-alpha(p,0))*(neu_alpha(p,0)-alpha(p,0));
        }

        double d_const = fabs(d(0,0));
        if(fabs(chi2-chisqrd)<1e-3 && d_const<10 && sqrt(sum0)<1e-3){
            fIteration = q;
            fConverged = true;
            break;
        }
        */
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
    fPull.ResizeTo(fN * cov_dim, fN * cov_dim);
    for (uint b = 0; b < (fN * cov_dim); b++)
        fPull(b, b) = -10000;

    if (true)
    {
        for (uint b = 0; b < (fN * cov_dim); b++)
        {
            double num = A0(b, 0) - alpha(b, 0);
            double dem = V0(b, b) - V(b, b);
            if (dem > 0) { fPull(b, b) = num / std::sqrt(dem); }
        }
    }

    updateDaughters();

    // return fConverged; // for number of iterations greater than 1
    return true; // for number of iterations equal to 1
}

HRefitCand HVertexFitter::getDaughter(int val)
{
    return fCands[val];
}

void HVertexFitter::updateDaughters()
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

void HVertexFitter::update()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand& cand = fCands[val];
        cand.update();
    }
}
