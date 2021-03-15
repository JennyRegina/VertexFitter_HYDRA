#include "hvertexfitter.h"

const size_t cov_dim = 5;

HVertexFitter::HVertexFitter(const std::vector<HRefitCand> &cands) : fCands(cands), fVerbose(0), fLearningRate(0.5), fNumIterations(5), fPrimaryVertexFound(false)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();

    y.ResizeTo(fN * cov_dim, 1);
    V.ResizeTo(fN * cov_dim, fN * cov_dim);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;

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

TVector3 HVertexFitter::findVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFitter::findVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double param_p1, param_theta1, param_phi1, param_R1, param_Z1;
    double param_p2, param_theta2, param_phi2, param_R2, param_Z2;

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
    energy_cand1 = sqrt(param_p1*param_p1+983.272*983.272);
    energy_cand2 = sqrt(param_p2*param_p2+139.570*139.570);

    double momentumAfterDecay = sqrt(energy_cand1*energy_cand1+2*energy_cand1*energy_cand2+energy_cand2*energy_cand2-1115.683*1115.683);

    //double momentumAfterDecay= param_p1 + param_p2;

    // Set the original phi angles
    //fPhi1Original=param_phi1;
    //fPhi2Original=param_phi2;

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

    // Create a Lambda candidate object.
    // As a first approximation it originates from x,y,z = 0,0,0

   std::cout << "Vertex: theta: " << fVertex.Theta() << " and phi: " << fVertex.Phi() << std::endl;

   if(fPrimaryVertexFound==false){
   	setLambdaCandidate(momentumAfterDecay, fVertex.Theta(), fVertex.Phi(), 0.0, 0.0, fVertex);
   }
   if(fPrimaryVertexFound==true){
	setLambdaCandidateFromPrimaryVtxInfo(momentumAfterDecay, fVertex.Theta(), fVertex.Phi(), 0.0, 0.0, fPrimaryVertex);
   }

    return fVertex;
}

std::vector<HRefitCand> HVertexFitter::UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &vertexPos)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFitter::UpdateTrackParameters() -----------" << std::endl;
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

    /*double R_secondary1=;
    double R_secondary2=;

    double Z_secondary1=;
    double Z_secondary2=; */

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

TMatrixD HVertexFitter::f_eval(const TMatrixD &m_iter)
{
    TMatrixD d;

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

    return d;
}

void HVertexFitter::setLambdaCandidate(double valMomentum, double valTheta, double valPhi, double valR, double valZ, TVector3 decayVertex){

//TODO make sure that all properties of the virtual can is set properly

fLambdaCandidate.setMomentum(valMomentum);
fLambdaCandidate.setTheta(TMath::RadToDeg() * valTheta);
fLambdaCandidate.setPhi(TMath::RadToDeg() * valPhi);
fLambdaCandidate.SetTheta(valTheta);
fLambdaCandidate.SetPhi(valPhi);
fLambdaCandidate.setR(valR);
fLambdaCandidate.setZ(valZ);

std::cout << "setLambdaCandidate, fLambdaCandidate: theta= " << fLambdaCandidate.getTheta() << " and phi = " << fLambdaCandidate.getPhi() << std::endl; 

// Calculate the covariance matrix for the Lambda Candidate

double x_vertex=decayVertex.X();
double y_vertex=decayVertex.Y();
double z_vertex=decayVertex.Z();


double sigma_x=33.39; // In mm
double sigma_y=26.70; // In mm
double sigma_z=44.92; // In mm

// Use coordinate transformation cartesian->polar to estimate error in theta and phi

// Calculate the error in theta
double r=sqrt(x_vertex*x_vertex+y_vertex*y_vertex+z_vertex*z_vertex);

double dtheta_dx=x_vertex*z_vertex/(r*r*r*sqrt(1-z_vertex/(r*r)));
double dtheta_dy=y_vertex*z_vertex/(r*r*r*sqrt(1-z_vertex/(r*r)));
double dtheta_dz=(1/r-z_vertex*z_vertex/(r*r*r))/sqrt(1-z_vertex*z_vertex/(r*r));

double sigma_theta=sqrt(dtheta_dx*dtheta_dx*sigma_x*sigma_x+dtheta_dy*dtheta_dy*sigma_y*sigma_y+dtheta_dz*dtheta_dz*sigma_z*sigma_z);

// Calculate the error in phi
double r_2D = sqrt(x_vertex*x_vertex+y_vertex*y_vertex);

double dphi_dx = -x_vertex*y_vertex/(sqrt(x_vertex*x_vertex/(r_2D*r_2D))*r_2D*r_2D*r_2D);
double dphi_dy = sqrt(x_vertex*x_vertex/(r_2D*r_2D))/r_2D;
// dphi_dz=0;

double sigma_phi = sqrt(dphi_dx*dphi_dx*sigma_x*sigma_x+dphi_dy*dphi_dy*sigma_y*sigma_y);

// Calculate the error in R
double dR_dx=x_vertex/r_2D;
double dR_dy=y_vertex/r_2D;

double sigma_R=sqrt(dR_dx*dR_dx*sigma_x*sigma_x+dR_dy*dR_dy*sigma_y*sigma_y);

fCovarianceLambda.ResizeTo(5,5);
fCovarianceLambda(0,0)=9999999;
fCovarianceLambda(1,1)=sigma_theta*sigma_theta;
fCovarianceLambda(2,2)=sigma_phi*sigma_phi;
fCovarianceLambda(3,3)=sigma_R*sigma_R;
fCovarianceLambda(4,4)=sigma_z*sigma_z;

}

TMatrixD HVertexFitter::Feta_eval(const TMatrixD &m_iter)
{

    TMatrixD H;

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

    return H;
}

bool HVertexFitter::fit()
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFitter::fit() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double lr = fLearningRate;
    TMatrixD alpha0(fN * cov_dim, 1), alpha(fN * cov_dim, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    TMatrixD alpha_original = alpha0;

    // Calculating the original covariance matrix that is not changed in the iterations
    TMatrixD V0_inv(V); // J.R New
    V0_inv.Invert();    // J.R New

    double chi2 = 1e6;
    TMatrixD D = Feta_eval(alpha);
    TMatrixD d = f_eval(alpha);

    for (int q = 0; q < fNumIterations; q++)
    {
        TMatrixD DT(D.GetNcols(), D.GetNrows());
        DT.Transpose(D);
        //TMatrixD VD = D * V * DT; Old
        TMatrixD VD = D * V0 * DT;
        VD.Invert();

        // TMatrixD delta_alpha = alpha - alpha0;
        TMatrixD delta_alpha = alpha_original - alpha;

        TMatrixD lambda = VD * D * delta_alpha + VD * d;
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fN * cov_dim, 1);
        // neu_alpha = alpha_original - lr * V * DT * lambda;
        neu_alpha = alpha_original - lr * V0 * DT * lambda; // J.R. New

        //double chisqrd = 0.;

        //Calculate new chi2
        TMatrixD chisqrd(1, 1);
        TMatrixD delta_alphaT(delta_alpha.GetNcols(), delta_alpha.GetNrows());
        delta_alphaT.Transpose(delta_alpha);
        TMatrixD two(1, 1);
        two(0, 0) = 2;
        chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * f_eval(neu_alpha);
        //chisqrd = f_eval(neu_alpha);
        //std::cout << "chisqrd = " << chisqrd << ", and chi2 = " << chi2 << std::endl;

        //chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * d;

        //for (int p = 0; p < lambda.GetNrows(); p++)
        //{
        //    chisqrd = lambdaT(0, p) * d(p, 0);
        //}

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
        //chi2 = chisqrd;

        /* if (fabs(chi2 - chisqrd(0, 0)) < 1)
        {

            fIteration = q;
            if (fVerbose > 0)
            {
                std::cout << "Iteration: " << q << ", Condition fullfilled!" << std::endl;
            }
            fConverged = true;
            chi2 = chisqrd(0, 0);
            alpha = neu_alpha;
            V = V0 - lr * V0 * DT * VD * D * V0;
            break;

        } */

        if (fVerbose > 0)
        {
            std::cout << "Iteration: " << q << std::endl;
            std::cout << "Printing d: " << std::endl; 
            d.Print();
        }

        fIteration = q;
        chi2 = chisqrd(0, 0);
        alpha0 = alpha;
        alpha = neu_alpha;
        //V = V - lr * V * DT * VD * D * V;
        V = V0 - lr * V0 * DT * VD * D * V0;
        D = Feta_eval(alpha);
        d = f_eval(alpha);


    }

    y = alpha;
    fChi2 = chi2;
    fNdf = 1; // The vertex constraint is a 1C fit
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
            if (dem > 0)
            {
                fPull(b, b) = num / std::sqrt(dem);
            }
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

void HVertexFitter::update()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand &cand = fCands[val];
        cand.update();
    }
}
