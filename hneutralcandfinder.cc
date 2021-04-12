#include "hneutralcandfinder.h"

HNeutralCandFinder::HNeutralCandFinder(const std::vector<HRefitCand> &cands) : fCands(cands), fVerbose(0), fPrimaryVertexFound(false), fUsePrimaryVertexInNeutralCandidateCalculation(false)
{
    double param_p1, param_p2;
    
    HRefitCand cand1 = cands[0];

    param_p1 = cand1.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate

    HRefitCand cand2 = cands[1];

    param_p2 = cand2.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate
    
    double energy_cand1, energy_cand2;
    energy_cand1 = sqrt(param_p1 * param_p1 + 983.272 * 983.272);
    energy_cand2 = sqrt(param_p2 * param_p2 + 139.570 * 139.570);

    double fMomentumAfterDecay = sqrt(energy_cand1 * energy_cand1 + 2 * energy_cand1 * energy_cand2 + energy_cand2 * energy_cand2 - 1115.683 * 1115.683);
    // If the primary vertex was found, the neutral mother candidate is created from the primary and decay vertex info

    /* if (fPrimaryVertexFound == false || fUsePrimaryVertexInNeutralCandidateCalculation == false)
    {
        setNeutralMotherCand(momentumAfterDecay, fVertex.Theta(), fVertex.Phi(), 0.0, 0.0, fVertex);
    }
    if (fPrimaryVertexFound == true && fUsePrimaryVertexInNeutralCandidateCalculation == true)
    {
        setNeutralMotherCandFromPrimaryVtxInfo(momentumAfterDecay, fPrimaryVertex, fVertex);
    } */
}

void HNeutralCandFinder::setNeutralMotherCand(double valMomentum, double valTheta, double valPhi, double valR, double valZ, TVector3 decayVertex)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HNeutralCandFinder::setNeutralMotherCandidate() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    //TODO make sure that all properties of the virtual cand is set properly
    //TODO use matrix notation to include the correlations in the errors

    fNeutralMotherCandidate.setMomentum(valMomentum);
    fNeutralMotherCandidate.setTheta(TMath::RadToDeg() * valTheta);
    fNeutralMotherCandidate.setPhi(TMath::RadToDeg() * valPhi);
    fNeutralMotherCandidate.SetTheta(valTheta);
    fNeutralMotherCandidate.SetPhi(valPhi);
    fNeutralMotherCandidate.setR(valR);
    fNeutralMotherCandidate.setZ(valZ);

    if (fVerbose > 0)
    {
        std::cout << "setNeutralMotherCandidate, fNeutralMotherCandidate: theta= " << fNeutralMotherCandidate.getTheta() << " and phi = " << fNeutralMotherCandidate.getPhi() << std::endl;
    }

    // Calculate the covariance matrix for the Lambda Candidate

    double x_vertex = decayVertex.X();
    double y_vertex = decayVertex.Y();
    double z_vertex = decayVertex.Z();

    // the errors below are estimated from difference distributions between reconstructed - MC truth for the vertex

    double sigma_x = 33.29; //old value: 33.39; // In mm
    double sigma_y = 27.09; //old value: 26.70; // In mm
    double sigma_z = 44.84; //old value: 44.92; // In mm

    // Use coordinate transformation cartesian->polar to estimate error in theta and phi

    // Calculate the error in theta
    double r = sqrt(x_vertex * x_vertex + y_vertex * y_vertex + z_vertex * z_vertex);

    double dtheta_dx = x_vertex * z_vertex / (r * r * r * sqrt(1 - z_vertex / (r * r)));
    double dtheta_dy = y_vertex * z_vertex / (r * r * r * sqrt(1 - z_vertex / (r * r)));
    double dtheta_dz = (1 / r - z_vertex * z_vertex / (r * r * r)) / sqrt(1 - z_vertex * z_vertex / (r * r));

    double sigma_theta = sqrt(dtheta_dx * dtheta_dx * sigma_x * sigma_x + dtheta_dy * dtheta_dy * sigma_y * sigma_y + dtheta_dz * dtheta_dz * sigma_z * sigma_z);

    // Calculate the error in phi
    double r_2D = sqrt(x_vertex * x_vertex + y_vertex * y_vertex);

    double dphi_dx = -x_vertex * y_vertex / (sqrt(x_vertex * x_vertex / (r_2D * r_2D)) * r_2D * r_2D * r_2D);
    double dphi_dy = sqrt(x_vertex * x_vertex / (r_2D * r_2D)) / r_2D;
    // dphi_dz=0;

    double sigma_phi = sqrt(dphi_dx * dphi_dx * sigma_x * sigma_x + dphi_dy * dphi_dy * sigma_y * sigma_y);

    // Calculate the error in R
    double dR_dx = x_vertex / r_2D;
    double dR_dy = y_vertex / r_2D;

    double sigma_R = sqrt(dR_dx * dR_dx * sigma_x * sigma_x + dR_dy * dR_dy * sigma_y * sigma_y);

    fCovarianceNeutralMother.ResizeTo(5, 5);
    fCovarianceNeutralMother(0, 0) = 9999999;
    fCovarianceNeutralMother(1, 1) = sigma_theta * sigma_theta;
    fCovarianceNeutralMother(2, 2) = sigma_phi * sigma_phi;
    fCovarianceNeutralMother(3, 3) = sigma_R * sigma_R;
    fCovarianceNeutralMother(4, 4) = sigma_z * sigma_z;
}

void HNeutralCandFinder::setNeutralMotherCandFromPrimaryVtxInfo(TVector3 primaryVertex, TVector3 decayVertex)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HNeutralCandFinder::setNeutralMotherCandFromPrimaryVtxInfo() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    // TODO: Set the other track parameters and the covariance matrix


    TVector3 vecPrimToDecayVertex = decayVertex - primaryVertex;

    double thetaPrimaryToSecondaryVertex, phiPrimaryToSecondaryVertex;

    thetaPrimaryToSecondaryVertex = vecPrimToDecayVertex.Theta();
    phiPrimaryToSecondaryVertex = vecPrimToDecayVertex.Phi();

    fNeutralMotherCandidate.setTheta(TMath::RadToDeg() * thetaPrimaryToSecondaryVertex);
    fNeutralMotherCandidate.setPhi(TMath::RadToDeg() * phiPrimaryToSecondaryVertex);
}
