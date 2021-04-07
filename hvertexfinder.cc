#include "hvertexfinder.h"

HVertexFinder::HVertexFinder(const std::vector<HRefitCand> &cands) : fCands(cands), fVerbose(0), fPrimaryVertexFound(false), fUsePrimaryVertexInNeutralCandidateCalculation(false)
{
}

TVector3 HVertexFinder::findVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::findVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double param_p1, param_theta1, param_phi1, param_R1, param_Z1;
    double param_p2, param_theta2, param_phi2, param_R2, param_Z2;

    if (cands.size() < 2)
    {

        std::cout << "WARNING: No vertex can be found, not enough charged particles in the event" << std::endl;
        fVertex.SetXYZ(-2000., -2000., -2000.);
        return fVertex;
    }

    HRefitCand cand1 = cands[0];

    param_p1 = cand1.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate
    param_theta1 = cand1.Theta();
    param_phi1 = cand1.Phi();
    param_R1 = cand1.getR();
    param_Z1 = cand1.getZ();

    HRefitCand cand2 = cands[1];

    param_p2 = cand2.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate
    param_theta2 = cand2.Theta();
    param_phi2 = cand2.Phi();
    param_R2 = cand2.getR();
    param_Z2 = cand2.getZ();

    double energy_cand1, energy_cand2;
    energy_cand1 = sqrt(param_p1 * param_p1 + 983.272 * 983.272);
    energy_cand2 = sqrt(param_p2 * param_p2 + 139.570 * 139.570);

    double momentumAfterDecay = sqrt(energy_cand1 * energy_cand1 + 2 * energy_cand1 * energy_cand2 + energy_cand2 * energy_cand2 - 1115.683 * 1115.683);

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    // Calculate the distance between the two tracks
    double dist = std::fabs((vtx_dir_1.Cross(vtx_dir_2)).Dot((vtx_base_1 - vtx_base_2)));
    fDistanceParticleToParticle = dist;

    if (fVerbose > 0)
    {
        std::cout << "Minimum distance between tracks: " << dist << std::endl;
        std::cout << " " << std::endl;
    }

    HGeomVector vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;

    // Direction vectors
    vtx_geom_dir_1.setX(std::sin(param_theta1) * std::cos(param_theta1));
    vtx_geom_dir_1.setY(std::sin(param_theta1) * std::sin(param_phi1));
    vtx_geom_dir_1.setZ(std::cos(param_theta1));
    vtx_geom_dir_2.setX(std::sin(param_theta2) * std::cos(param_theta2));
    vtx_geom_dir_2.setY(std::sin(param_theta2) * std::sin(param_phi2));
    vtx_geom_dir_2.setZ(std::cos(param_theta2));

    // Base vectors
    vtx_geom_base_1.setX(param_R1 * std::cos(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.setY(param_R1 * std::sin(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.setZ(param_Z1);
    vtx_geom_base_2.setX(param_R2 * std::cos(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.setY(param_R2 * std::sin(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.setZ(param_Z2);

    if (fVerbose > 0)
    {
        std::cout << " " << std::endl;
        std::cout << " ---------------------------------------" << std::endl;
        std::cout << " " << std::endl;
        std::cout << "Printing base vector 1: " << std::endl;
        vtx_geom_base_1.print();

        std::cout << " " << std::endl;
        std::cout << "Printing direction vector 1: " << std::endl;
        vtx_geom_dir_1.print();
        std::cout << " " << std::endl;

        std::cout << " " << std::endl;
        std::cout << "Printing base vector 2: " << std::endl;
        vtx_geom_base_2.print();

        std::cout << " " << std::endl;
        std::cout << "Printing direction vector 2: " << std::endl;
        vtx_geom_dir_2.print();
        std::cout << " " << std::endl;
        std::cout << " ---------------------------------------" << std::endl;
        std::cout << " " << std::endl;
    }

    HGeomVector vertex = HParticleTool::calculatePointOfClosestApproach(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_dir_2);

    fVertex.SetXYZ(vertex.X(), vertex.Y(), vertex.Z());

    double distanceFromParticleToVertex_1 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_1, vtx_geom_dir_1, vertex);
    double distanceFromParticleToVertex_2 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_2, vtx_geom_dir_2, vertex);

    HGeomVector originVertex;
    originVertex.setX(0.0);
    originVertex.setY(0.0);
    originVertex.setZ(0.0);

    double distanceFromParticleToOrigin_1 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_1, vtx_geom_dir_1, originVertex);
    double distanceFromParticleToOrigin_2 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_2, vtx_geom_dir_2, originVertex);

    fDistParticle1Vertex = distanceFromParticleToVertex_1;
    fDistParticle2Vertex = distanceFromParticleToVertex_2;
    fDistParticle1Origin = distanceFromParticleToOrigin_1;
    fDistParticle2Origin = distanceFromParticleToOrigin_2;

    // Find the primary vertex

    if (fVerbose > 0)
    {
        std::cout << "Vertex: theta: " << fVertex.Theta() << " and phi: " << fVertex.Phi() << std::endl;
    }

    fPrimaryVertexFound = false;

    if (cands.size() > 2)
    {
        findPrimaryVertex(cands);
    }

    if (fPrimaryVertexFound == true)
    {
        calculateVertexProperties(fPrimaryVertex, fVertex);
    }

    // If the primary vertex was found, the neutral mother candidate is created from the primary and decay vertex info

    if (fPrimaryVertexFound == false || fUsePrimaryVertexInNeutralCandidateCalculation == false)
    {
        setNeutralMotherCandidate(momentumAfterDecay, fVertex.Theta(), fVertex.Phi(), 0.0, 0.0, fVertex);
    }
    if (fPrimaryVertexFound == true && fUsePrimaryVertexInNeutralCandidateCalculation == true)
    {
        setNeutralMotherCandidateFromPrimaryVtxInfo(momentumAfterDecay, fPrimaryVertex, fVertex);
    }

    return fVertex;
}

std::vector<HRefitCand> HVertexFinder::UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &vertexPos)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::UpdateTrackParameters() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double param_theta1, param_phi1, param_R1, param_Z1;
    double param_theta2, param_phi2, param_R2, param_Z2;

    HRefitCand cand1 = cands[0];

    //param_p_inv1 = 1. / cand1.P();
    param_theta1 = cand1.Theta();
    param_phi1 = cand1.Phi();
    param_R1 = cand1.getR();
    param_Z1 = cand1.getZ();

    HRefitCand cand2 = cands[1];

    //param_p_inv2 = 1. / cand2.P();
    param_theta2 = cand2.Theta();
    param_phi2 = cand2.Phi();
    param_R2 = cand2.getR();
    param_Z2 = cand2.getZ();

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    //Vectors pointing from vertex to POCA to Beam axis
    TVector3 vtx_base_1_updated, vtx_base_2_updated;
    vtx_base_1_updated = vtx_base_1 - vertexPos;
    vtx_base_2_updated = vtx_base_2 - vertexPos;

    if (fVerbose > 0)
    {
        std::cout << " " << std::endl;
        std::cout << "Position of base vector 1: x=" << vtx_base_1.X() << ", y=" << vtx_base_1.Y() << ", z=" << vtx_base_1.Z() << std::endl;
        std::cout << "Position of updated base vector 1: x=" << vtx_base_1_updated.X() << ", y=" << vtx_base_1_updated.Y() << ", z=" << vtx_base_1_updated.Z() << std::endl;
        std::cout << "Theta, original base vector 1: " << vtx_base_1.Theta() << ", Phi, original base vector: " << vtx_base_1.Phi() << std::endl;
        std::cout << "Theta, updated base vector 1: " << vtx_base_1_updated.Theta() << ", Phi, updated base vector: " << vtx_base_1_updated.Phi() << std::endl;
        std::cout << " " << std::endl;

        std::cout << " " << std::endl;
        std::cout << "Position of base vector 2: x=" << vtx_base_2.X() << ", y=" << vtx_base_2.Y() << ", z=" << vtx_base_2.Z() << std::endl;
        std::cout << "Position of updated base vector 2: x=" << vtx_base_2_updated.X() << ", y=" << vtx_base_2_updated.Y() << ", z=" << vtx_base_2_updated.Z() << std::endl;
        std::cout << "Theta, original base vector 2: " << vtx_base_2.Theta() << ", Phi, original base vector: " << vtx_base_2.Phi() << std::endl;
        std::cout << "Theta, updated base vector 2: " << vtx_base_2_updated.Theta() << ", Phi, updated base vector: " << vtx_base_2_updated.Phi() << std::endl;
        std::cout << " " << std::endl;
    }

    double theta_secondary1 = vtx_base_1_updated.Theta();
    double theta_secondary2 = vtx_base_2_updated.Theta();

    double phi_secondary1 = vtx_base_1_updated.Phi();
    double phi_secondary2 = vtx_base_2_updated.Phi();

    TVector3 vtx_dir_1_updated, vtx_dir_2_updated;

    vtx_dir_1_updated.SetXYZ(std::sin(theta_secondary1) * std::cos(phi_secondary1),
                             std::sin(theta_secondary1) * std::sin(phi_secondary1),
                             std::cos(theta_secondary1));
    vtx_dir_2_updated.SetXYZ(std::sin(theta_secondary2) * std::cos(phi_secondary2),
                             std::sin(theta_secondary2) * std::sin(phi_secondary2),
                             std::cos(theta_secondary2));

    // Calculate the distance between the two tracks
    double dist_new = std::fabs((vtx_dir_1_updated.Cross(vtx_dir_2_updated)).Dot((vtx_base_1_updated - vtx_base_2_updated)));

    if (fVerbose > 0)
    {
        std::cout << "Minimum NEW distance between tracks: " << dist_new << std::endl;
        std::cout << " " << std::endl;
    }

    if (fVerbose > 0)
    {
        std::cout << "Before update " << std::endl;

        std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
        std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
        std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
        std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
    }

    cand1.SetTheta(theta_secondary1);
    cand1.SetPhi(phi_secondary1);

    cand2.SetTheta(theta_secondary2);
    cand2.SetPhi(phi_secondary2);

    if (fVerbose > 0)
    {
        std::cout << "After update " << std::endl;

        std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
        std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
        std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
        std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
    }

    std::vector<HRefitCand> newCands;
    newCands.clear();
    newCands.push_back(cand1);
    newCands.push_back(cand2);

    return newCands;
}

TVector3 HVertexFinder::findPrimaryVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::findPrimaryVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double param_theta1, param_phi1, param_R1, param_Z1;
    double param_theta2, param_phi2, param_R2, param_Z2;

    // All found protons in the event
    HRefitCand primaryCand1 = cands[0];

    param_theta1 = primaryCand1.Theta();
    param_phi1 = primaryCand1.Phi();
    param_R1 = primaryCand1.getR();
    param_Z1 = primaryCand1.getZ();

    HRefitCand primaryCand2 = cands[2];
    param_theta2 = primaryCand2.Theta();
    param_phi2 = primaryCand2.Phi();
    param_R2 = primaryCand2.getR();
    param_Z2 = primaryCand2.getZ();

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    HGeomVector vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;

    // Direction vectors
    vtx_geom_dir_1.setX(std::sin(param_theta1) * std::cos(param_theta1));
    vtx_geom_dir_1.setY(std::sin(param_theta1) * std::sin(param_phi1));
    vtx_geom_dir_1.setZ(std::cos(param_theta1));
    vtx_geom_dir_2.setX(std::sin(param_theta2) * std::cos(param_theta2));
    vtx_geom_dir_2.setY(std::sin(param_theta2) * std::sin(param_phi2));
    vtx_geom_dir_2.setZ(std::cos(param_theta2));

    // Base vectors
    vtx_geom_base_1.setX(param_R1 * std::cos(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.setY(param_R1 * std::sin(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.setZ(param_Z1);
    vtx_geom_base_2.setX(param_R2 * std::cos(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.setY(param_R2 * std::sin(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.setZ(param_Z2);

    HGeomVector primaryVertex = HParticleTool::calculatePointOfClosestApproach(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_dir_2);

    fPrimaryVertex.SetXYZ(primaryVertex.X(), primaryVertex.Y(), primaryVertex.Z());

    // if the primary vertex was not found, each coordinate x,y,z is set to -20000
    if (primaryVertex.X() != -2000)
    {
        fPrimaryVertexFound = true;
    }

    return fPrimaryVertex;
}

void HVertexFinder::setNeutralMotherCandidate(double valMomentum, double valTheta, double valPhi, double valR, double valZ, TVector3 decayVertex)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::setNeutralMotherCandidate() -----------" << std::endl;
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

    double sigma_x = 33.39; // In mm
    double sigma_y = 26.70; // In mm
    double sigma_z = 44.92; // In mm

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

void HVertexFinder::setNeutralMotherCandidateFromPrimaryVtxInfo(double valMomentum, TVector3 primaryVertex, TVector3 decayVertex)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::setNeutralMotherCandidateFromPrimaryVtxInfo() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    // TODO: Set the other track parameters and the covariance matrix

    calculateVertexProperties(primaryVertex, decayVertex);

    double thetaPrimaryToSecondaryVertex, phiPrimaryToSecondaryVertex;

    thetaPrimaryToSecondaryVertex = fVecPrimToDecayVertex.Theta();
    phiPrimaryToSecondaryVertex = fVecPrimToDecayVertex.Phi();

    fNeutralMotherCandidate.setTheta(TMath::RadToDeg() * thetaPrimaryToSecondaryVertex);
    fNeutralMotherCandidate.setPhi(TMath::RadToDeg() * phiPrimaryToSecondaryVertex);
}

void HVertexFinder::calculateVertexProperties(TVector3 primaryVertex, TVector3 decayVertex)
{

    fVecPrimToDecayVertex = decayVertex - primaryVertex;

    fDistPrimToDecayVertex = sqrt((decayVertex.X() - primaryVertex.X()) * (decayVertex.X() - primaryVertex.X()) + (decayVertex.Y() - primaryVertex.Y()) * (decayVertex.Y() - primaryVertex.Y()) + (decayVertex.Z() - primaryVertex.Z()) * (decayVertex.Z() - primaryVertex.Z()));

    if (decayVertex.Z() > primaryVertex.Z())
    {
        fPrimaryVertexIsBetforeDecayVertex = true;
    }
    else
    {
        fPrimaryVertexIsBetforeDecayVertex = false;
    }

    double R_primaryVertex, R_decayVertex;

    R_primaryVertex = sqrt(primaryVertex.X() * primaryVertex.X() + primaryVertex.Y() * primaryVertex.Y());
    R_decayVertex = sqrt(decayVertex.X() * decayVertex.X() + decayVertex.Y() * decayVertex.Y());

    if (R_primaryVertex < R_decayVertex)
    {
        fPrimaryVertexIsInsideDecayVertex = true;
    }
    else
    {
        fPrimaryVertexIsInsideDecayVertex = false;
    }
}
