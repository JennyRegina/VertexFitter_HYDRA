#include "hades.h"
#include "hcategorymanager.h"
#include "henergylosscorrpar.h"
#include "hhistmap.h"
#include "hloop.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"

#include "hcategory.h"
#include "hlinearcategory.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlegeantpair.h"
#include "hparticlepair.h"
#include "hrichhit.h"
#include "hrichhitsim.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"

#include "hparticlecutrange.h"
#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"

#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "hvertexfitter.h"
#include "hvertexfinder.h"

using namespace std;
using namespace Particle;

void FillData(HParticleCand *cand, HRefitCand &outcand, double arr[],
              double mass)
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand.SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand.setR(cand->getR());
    outcand.setZ(cand->getZ());
    outcand.setCovariance(cov);
}

Bool_t selectHadrons(HParticleCand *pcand)
{
    // build in selection function for hadron candidates.
    // Requires besides an RK + META and fitted
    // inner+outer segment.

    Bool_t test = kFALSE;
    if (pcand->isFlagAND(4, Particle::kIsAcceptedHitInnerMDC,
                         Particle::kIsAcceptedHitOuterMDC,
                         Particle::kIsAcceptedHitMETA,
                         Particle::kIsAcceptedRK) &&
        pcand->getInnerSegmentChi2() > 0 && pcand->getChi2() < 10000 // RK
    )
        test = kTRUE;

    if (!test)
        return kFALSE;

    if (test)
        test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE;

    return test;
}

Int_t analysisVertexFinder(TString infileList = "pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents = 100000)
{

    // Min input: pp_pKlambda_100000evts1_dst_apr12.root

    TStopwatch timer;
    timer.Start();

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile *outfile = new TFile("vertexfit_WithoutCombBackg.root", "recreate");

    TH1F *h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F *h02 = new TH1F("hChi2", "", 100, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" Counts ");

    TH1F *h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" Counts ");

    TH1F *h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" Events ");

    //TH1F *h17 = new TH1F("hLambdaMassPostFitCut", "", 100, 1070, 1170);
    //h17->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    //h17->SetYTitle(" Events ");

    TH1F *h05 = new TH1F("hPullPInv", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" Counts ");

    TH1F *h06 = new TH1F("hPullTheta", "", 100, -5, 5);
    h06->SetXTitle("Pull(#theta)");
    h06->SetYTitle(" Counts ");

    TH1F *h07 = new TH1F("hPullPhi", "", 100, -5, 5);
    h07->SetXTitle("Pull(#phi)");
    h07->SetYTitle(" Counts ");

    TH1F *h08 = new TH1F("hPullR", "", 100, -5, 5);
    h08->SetXTitle("Pull(R)");
    h08->SetYTitle(" Counts ");

    TH1F *h09 = new TH1F("hPullZ", "", 100, -5, 5);
    h09->SetXTitle("Pull(Z)");
    h09->SetYTitle(" Counts ");

    /*  TH1F *h10 = new TH1F("hSuccessfulConvergence", "", 100, 1070, 1170);
    h10->SetXTitle("M_{p#pi^{-}} [MeV/c^{2}]");
    h10->SetYTitle(" Counts "); */

    TH1F *h10 = new TH1F("hSuccessfulConvergence", "", 2, 0, 2);
    h10->SetXTitle("Fit has converged");
    h10->SetYTitle(" Counts ");

    TH1F *h11 = new TH1F("hIterations", "", 100, 0, 20);
    h11->SetXTitle("Number of Iterations");
    h11->SetYTitle(" Counts ");

    // ------------------------- Histos after probability cut    --------------------

    TH1F *h12 = new TH1F("hPullPInvCut", "", 100, -5, 5);
    h12->SetXTitle("Pull(1/P_{p})");
    h12->SetYTitle(" Counts ");

    TH1F *h13 = new TH1F("hPullThetaCut", "", 100, -5, 5);
    h13->SetXTitle("Pull(#theta)");
    h13->SetYTitle(" Counts ");

    TH1F *h14 = new TH1F("hPullPhiCut", "", 100, -5, 5);
    h14->SetXTitle("Pull(#phi)");
    h14->SetYTitle(" Counts ");

    TH1F *h15 = new TH1F("hPullRCut", "", 100, -5, 5);
    h15->SetXTitle("Pull(R)");
    h15->SetYTitle(" Counts ");

    TH1F *h16 = new TH1F("hPullZCut", "", 100, -5, 5);
    h16->SetXTitle("Pull(Z)");
    h16->SetYTitle(" Counts ");

    // -------------------- Pulls Pions -----------------------

    TH1F *h17 = new TH1F("hPullThetaPion", "", 100, -5, 5);
    h17->SetXTitle("Pull(#theta)");
    h17->SetYTitle(" Counts ");

    TH1F *h18 = new TH1F("hPullPhiPion", "", 100, -5, 5);
    h18->SetXTitle("Pull(#phi)");
    h18->SetYTitle(" Counts ");

    TH1F *h19 = new TH1F("hPullRPion", "", 100, -5, 5);
    h19->SetXTitle("Pull(R)");
    h19->SetYTitle(" Counts ");

    TH1F *h20 = new TH1F("hPullZPion", "", 100, -5, 5);
    h20->SetXTitle("Pull(Z)");
    h20->SetYTitle(" Counts ");

    // ------------------------- Histos after probability cut    --------------------

    TH1F *h21 = new TH1F("hPullThetaCutPion", "", 100, -5, 5);
    h21->SetXTitle("Pull(#theta)");
    h21->SetYTitle(" Counts ");

    TH1F *h22 = new TH1F("hPullPhiCutPion", "", 100, -5, 5);
    h22->SetXTitle("Pull(#phi)");
    h22->SetYTitle(" Counts ");

    TH1F *h23 = new TH1F("hPullRCutPion", "", 100, -5, 5);
    h23->SetXTitle("Pull(R)");
    h23->SetYTitle(" Counts ");

    TH1F *h24 = new TH1F("hPullZCutPion", "", 100, -5, 5);
    h24->SetXTitle("Pull(Z)");
    h24->SetYTitle(" Counts ");

    // -----------------------------------------------------------------------
    // ------------------------- Histos after convergence    --------------------

    TH1F *h25 = new TH1F("hPullPInvConverged", "", 100, -5, 5);
    h25->SetXTitle("Pull(1/P_{p})");
    h25->SetYTitle(" Counts ");

    TH1F *h26 = new TH1F("hPullThetaConverged", "", 100, -5, 5);
    h26->SetXTitle("Pull(#theta)");
    h26->SetYTitle(" Counts ");

    TH1F *h27 = new TH1F("hPullPhiConverged", "", 100, -5, 5);
    h27->SetXTitle("Pull(#phi)");
    h27->SetYTitle(" Counts ");

    TH1F *h28 = new TH1F("hPullRConverged", "", 100, -5, 5);
    h28->SetXTitle("Pull(R)");
    h28->SetYTitle(" Counts ");

    TH1F *h29 = new TH1F("hPullZConverged", "", 100, -5, 5);
    h29->SetXTitle("Pull(Z)");
    h29->SetYTitle(" Counts ");

    TH1F *h30 = new TH1F("hPullThetaPionConverged", "", 100, -5, 5);
    h30->SetXTitle("Pull(#theta)");
    h30->SetYTitle(" Counts ");

    TH1F *h31 = new TH1F("hPullPhiPioConvergedn", "", 100, -5, 5);
    h31->SetXTitle("Pull(#phi)");
    h31->SetYTitle(" Counts ");

    TH1F *h32 = new TH1F("hPullRPionConverged", "", 100, -5, 5);
    h32->SetXTitle("Pull(R)");
    h32->SetYTitle(" Counts ");

    TH1F *h33 = new TH1F("hPullZPionConverged", "", 100, -5, 5);
    h33->SetXTitle("Pull(Z)");
    h33->SetYTitle(" Counts ");

    // -----------------------------------------------------------------------
    // ------- Vertex histograms pre fit -----------------------

    TH1F *hVertexXPreFit = new TH1F("hVertexXPreFit", "", 1000, -100, 100);
    hVertexXPreFit->SetXTitle("Vertex, X / mm");
    hVertexXPreFit->SetYTitle(" Counts ");

    TH1F *hVertexYPreFit = new TH1F("hVertexYPreFit", "", 1000, -100, 100);
    hVertexYPreFit->SetXTitle("Vertex, Y / mm");
    hVertexYPreFit->SetYTitle(" Counts ");

    TH1F *hVertexZPreFit = new TH1F("hVertexZPreFit", "", 1000, -60, 1000);
    hVertexZPreFit->SetXTitle("Vertex, Z / mm");
    hVertexZPreFit->SetYTitle(" Counts ");

    TH2F *hVertexPreFit = new TH2F("hVertexPreFit", "", 1000, -100, 100, 1000, -100, 100);
    hVertexPreFit->SetXTitle("Vertex, Z / mm");
    hVertexPreFit->SetYTitle("Vertex, R / mm");

    TH1F *hPrimaryVertexXPreFit = new TH1F("hPrimaryVertexXPreFit", "", 1000, -100, 100);
    hPrimaryVertexXPreFit->SetXTitle("Vertex, X / mm");
    hPrimaryVertexXPreFit->SetYTitle(" Counts ");

    TH1F *hPrimaryVertexYPreFit = new TH1F("hPrimaryVertexYPreFit", "", 1000, -100, 100);
    hPrimaryVertexYPreFit->SetXTitle("Vertex, Y / mm");
    hPrimaryVertexYPreFit->SetYTitle(" Counts ");

    TH1F *hPrimaryVertexZPreFit = new TH1F("hPrimaryVertexZPreFit", "", 1000, -60, 1000);
    hPrimaryVertexZPreFit->SetXTitle("Vertex, Z / mm");
    hPrimaryVertexZPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceBetweenProtonAndPionPreFit = new TH1F("hDistanceBetweenProtonAndPionPreFit", "", 500, 0, 50);
    hDistanceBetweenProtonAndPionPreFit->SetXTitle("Distance between particles / mm");
    hDistanceBetweenProtonAndPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexProtonPreFit = new TH1F("hDistanceToVertexProtonPreFit", "", 1000, 0, 100);
    hDistanceToVertexProtonPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexProtonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexPionPreFit = new TH1F("hDistanceToVertexPionPreFit", "", 1000, 0, 100);
    hDistanceToVertexPionPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginProtonPreFit = new TH1F("hDistanceToOriginProtonPreFit", "", 1000, 0, 100);
    hDistanceToOriginProtonPreFit->SetXTitle("Distance to origin / mm");
    hDistanceToOriginProtonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginPionPreFit = new TH1F("hDistanceToOriginPionPreFit", "", 1000, 0, 100);
    hDistanceToOriginPionPreFit->SetXTitle("Distance to Origin / mm");
    hDistanceToOriginPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginSumPreFit = new TH1F("hDistanceToOriginSumPreFit", "", 1000, 0, 100);
    hDistanceToOriginSumPreFit->SetXTitle("Distance to Origin Sum / mm");
    hDistanceToOriginSumPreFit->SetYTitle(" Counts ");

    // --------- Vertex histograms post fit --------------------

    TH1F *hVertexXPostFit = new TH1F("hVertexXPostFit", "", 1000, -100, 100);
    hVertexXPostFit->SetXTitle("Vertex, X / mm");
    hVertexXPostFit->SetYTitle(" Counts ");

    TH1F *hVertexYPostFit = new TH1F("hVertexYPostFit", "", 1000, -100, 100);
    hVertexYPostFit->SetXTitle("Vertex, Y / mm");
    hVertexYPostFit->SetYTitle(" Counts ");

    TH1F *hVertexZPostFit = new TH1F("hVertexZPostFit", "", 1000, -100, 100);
    hVertexZPostFit->SetXTitle("Vertex, Z / mm");
    hVertexZPostFit->SetYTitle(" Counts ");

    TH2F *hVertexPostFit = new TH2F("hVertexPostFit", "", 1000, -100, 100, 1000, -100, 100);
    hVertexPostFit->SetXTitle("Vertex, Z / mm");
    hVertexPostFit->SetYTitle("Vertex, R /mm");

    TH1F *hDistanceBetweenProtonAndPionPostFit = new TH1F("hDistanceBetweenProtonAndPionPostFit", "", 1000, 0, 100);
    hDistanceBetweenProtonAndPionPostFit->SetXTitle("Distance between particles / mm");
    hDistanceBetweenProtonAndPionPostFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexProtonPostFit = new TH1F("hDistanceToVertexProtonPostFit", "", 1000, 0, 100);
    hDistanceToVertexProtonPostFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexProtonPostFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexPionPostFit = new TH1F("hDistanceToVertexPionPostFit", "", 1000, 0, 100);
    hDistanceToVertexPionPostFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexPionPostFit->SetYTitle(" Counts ");

    //---------------------- MC (Geant) information ------------------------

    TH1F *hGeantVertexXAll = new TH1F("hGeantVertexXAll", "", 1000, -100, 100);
    hGeantVertexXAll->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexYAll = new TH1F("hGeantVertexYAll", "", 1000, -100, 100);
    hGeantVertexYAll->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexZAll = new TH1F("hGeantVertexZAll", "", 1000, -100, 100);
    hGeantVertexZAll->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexXLambda = new TH1F("hGeantVertexXLambda", "", 1000, -100, 100);
    hGeantVertexXLambda->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexYLambda = new TH1F("hGeantVertexYLambda", "", 1000, -100, 100);
    hGeantVertexYLambda->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexZLambda = new TH1F("hGeantVertexZLambda", "", 1000, -100, 100);
    hGeantVertexZLambda->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand1 = new TH1F("hGeantVertexXCand1", "", 1000, -100, 100);
    hGeantVertexXCand1->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand1 = new TH1F("hGeantVertexYCand1", "", 1000, -100, 100);
    hGeantVertexYCand1->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand1 = new TH1F("hGeantVertexZCand1", "", 1000, -60, 1000);
    hGeantVertexZCand1->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand2 = new TH1F("hGeantVertexXCand2", "", 1000, -100, 100);
    hGeantVertexXCand2->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand2->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand2 = new TH1F("hGeantVertexYCand2", "", 1000, -100, 100);
    hGeantVertexYCand2->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand2->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand2 = new TH1F("hGeantVertexZCand2", "", 1000, -60, 1000);
    hGeantVertexZCand2->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand2->SetYTitle(" Counts ");

    TH2F *hGeantVertex = new TH2F("hGeantVertex", "", 1000, -100, 100, 1000, -100, 100);
    hGeantVertex->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertex->SetYTitle("Geant Vertex, R / mm");

    TH1F *hVertexXDiff = new TH1F("hVertexXDiff", "", 1000, -100, 100);
    hVertexXDiff->SetXTitle("Difference Vertex, X / mm");
    hVertexXDiff->SetYTitle(" Counts ");

    TH1F *hVertexYDiff = new TH1F("hVertexYDiff", "", 1000, -100, 100);
    hVertexYDiff->SetXTitle("Difference Vertex, Y / mm");
    hVertexYDiff->SetYTitle(" Counts ");

    TH1F *hVertexZDiff = new TH1F("hVertexZDiff", "", 1000, -100, 100);
    hVertexZDiff->SetXTitle("Difference Vertex, Z / mm");
    hVertexZDiff->SetYTitle(" Counts ");

    // -------------- Geant Hitograms --------------------

    // ---- Momentum -----
    TH1F *hGeantTotMomentumProtons = new TH1F("hGeantTotMomentumProtons", "", 1000, 0, 3000);
    hGeantTotMomentumProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumProtons = new TH1F("hGeantXMomentumProtons", "", 1000, -1000, 1000);
    hGeantXMomentumProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumProtons = new TH1F("hGeantYMomentumProtons", "", 1000, -1000, 1000);
    hGeantYMomentumProtons->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumProtons = new TH1F("hGeantZMomentumProtons", "", 1000, -100, 3000);
    hGeantZMomentumProtons->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumProtons->SetYTitle(" Counts ");

    TH1F *hGeantTotMomentumPions = new TH1F("hGeantTotMomentumPions", "", 1000, 0, 3000);
    hGeantTotMomentumPions->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumPions = new TH1F("hGeantXMomentumPions", "", 1000, -1000, 1000);
    hGeantXMomentumPions->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumPions = new TH1F("hGeantYMomentumPions", "", 1000, -1000, 1000);
    hGeantYMomentumPions->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumPions = new TH1F("hGeantZMomentumPions", "", 1000, -100, 3000);
    hGeantZMomentumPions->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumPions->SetYTitle(" Counts ");

    // ------- Reconstructed quantities -----------
    TH1F *hRecoRProtons = new TH1F("hRecoRProtons", "", 1000, -100, 100);
    hRecoRProtons->SetXTitle("R / mm");
    hRecoRProtons->SetYTitle(" Counts ");
    TH1F *hRecoZProtons = new TH1F("hRecoZProtons", "", 1000, -100, 100);
    hRecoZProtons->SetXTitle("Z / mm");
    hRecoZProtons->SetYTitle(" Counts ");
    TH1F *hRecoThetaProtons = new TH1F("hRecoThetaProtons", "", 500, 0, 3);
    hRecoThetaProtons->SetXTitle("#theta / rad");
    hRecoThetaProtons->SetYTitle(" Counts ");
    TH1F *hRecoPhiProtons = new TH1F("hRecoPhiProtons", "", 500, -4, 4);
    hRecoPhiProtons->SetXTitle("#phi / rad");
    hRecoPhiProtons->SetYTitle(" Counts ");
    TH1F *hRecoMomentumProtons = new TH1F("hRecoMomentumProtons", "", 1000, 0, 3000);
    hRecoMomentumProtons->SetXTitle("Momentum / MeV/c");
    hRecoMomentumProtons->SetYTitle(" Counts ");
    TH1F *hRecoBetaProtons = new TH1F("hRecoBetaProtons", "", 100, -0.1, 1);
    hRecoBetaProtons->SetXTitle("#beta");
    hRecoBetaProtons->SetYTitle(" Counts ");

    TH1F *hRecoRPions = new TH1F("hRecoRPions", "", 1000, -100, 100);
    hRecoRPions->SetXTitle("R / mm");
    hRecoRPions->SetYTitle(" Counts ");
    TH1F *hRecoZPions = new TH1F("hRecoZPions", "", 1000, -100, 100);
    hRecoZPions->SetXTitle("Z / mm");
    hRecoZPions->SetYTitle(" Counts ");
    TH1F *hRecoThetaPions = new TH1F("hRecoThetaPions", "", 500, 0, 3);
    hRecoThetaPions->SetXTitle("#theta / rad");
    hRecoThetaPions->SetYTitle(" Counts ");
    TH1F *hRecoPhiPions = new TH1F("hRecoPhiPions", "", 500, -4, 4);
    hRecoPhiPions->SetXTitle("#phi / rad");
    hRecoPhiPions->SetYTitle(" Counts ");
    TH1F *hRecoMomentumPions = new TH1F("hRecoMomentumPions", "", 1000, 0, 3000);
    hRecoMomentumPions->SetXTitle("Momentum / MeV/c");
    hRecoMomentumPions->SetYTitle(" Counts ");
    TH1F *hRecoBetaPions = new TH1F("hRecoBetaPions", "", 100, -0.1, 1);
    hRecoBetaPions->SetXTitle("#beta");
    hRecoBetaPions->SetYTitle(" Counts ");

    // ---------------- LAMBDA PLOTS --------------------
    TH1F *hMomLambda = new TH1F("hMomLambda", "", 1000, 0, 3000);
    hMomLambda->SetXTitle("Momentum / MeV/c");
    hMomLambda->SetYTitle(" Counts ");
    TH1F *hRecoThetaLambda = new TH1F("hRecoThetaLambda", "", 500, 0, 3);
    hRecoThetaLambda->SetXTitle("#theta / rad");
    hRecoThetaLambda->SetYTitle(" Counts ");
    TH1F *hRecoPhiLambda = new TH1F("hRecoPhiLambda", "", 500, -4, 4);
    hRecoPhiLambda->SetXTitle("#phi / rad");
    hRecoPhiLambda->SetYTitle(" Counts ");

    TH1F *hErrorRLambda = new TH1F("hErrorRLambda", "", 1000, -100, 100);
    hErrorRLambda->SetXTitle("R / mm");
    hErrorRLambda->SetYTitle(" Counts ");
    TH1F *hErrorZLambda = new TH1F("hErrorZLambda", "", 1000, -100, 100);
    hErrorZLambda->SetXTitle("Z / mm");
    hErrorZLambda->SetYTitle(" Counts ");
    TH1F *hErrorThetaLambda = new TH1F("hErrorThetaLambda", "", 500, 0, 3);
    hErrorThetaLambda->SetXTitle("#theta / rad");
    hErrorThetaLambda->SetYTitle(" Counts ");
    TH1F *hErrorPhiLambda = new TH1F("hErrorPhiLambda", "", 500, -4, 4);
    hErrorPhiLambda->SetXTitle("#phi / rad");
    hErrorPhiLambda->SetYTitle(" Counts ");

    // ---------------- LAMBDA PLOTS CUT --------------------
    TH1F *hMomLambdaCut = new TH1F("hMomLambdaCut", "", 1000, 0, 3000);
    hMomLambdaCut->SetXTitle("Momentum / MeV/c");
    hMomLambdaCut->SetYTitle(" Counts ");
    TH1F *hRecoThetaLambdaCut = new TH1F("hRecoThetaLambdaCut", "", 500, 0, 3);
    hRecoThetaLambdaCut->SetXTitle("#theta / rad");
    hRecoThetaLambdaCut->SetYTitle(" Counts ");
    TH1F *hRecoPhiLambdaCut = new TH1F("hRecoPhiLambdaCut", "", 500, -4, 4);
    hRecoPhiLambdaCut->SetXTitle("#phi / rad");
    hRecoPhiLambdaCut->SetYTitle(" Counts ");

    TH1F *hErrorRLambdaCut = new TH1F("hErrorRLambdaCut", "", 1000, -100, 100);
    hErrorRLambdaCut->SetXTitle("R / mm");
    hErrorRLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorZLambdaCut = new TH1F("hErrorZLambdaCut", "", 1000, -100, 100);
    hErrorZLambdaCut->SetXTitle("Z / mm");
    hErrorZLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorThetaLambdaCut = new TH1F("hErrorThetaLambdaCut", "", 500, 0, 3);
    hErrorThetaLambdaCut->SetXTitle("#theta / rad");
    hErrorThetaLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorPhiLambdaCut = new TH1F("hErrorPhiLambdaCut", "", 500, -4, 4);
    hErrorPhiLambdaCut->SetXTitle("#phi / rad");
    hErrorPhiLambdaCut->SetYTitle(" Counts ");

    // Vertex prim and dicay info
    TH1F *hDistPrimToDecayVertex = new TH1F("hDistPrimToDecayVertex", "", 500, 0, 100);
    hDistPrimToDecayVertex->SetXTitle(" Distance Between Primary and Decay Vertex ");
    hDistPrimToDecayVertex->SetYTitle(" Counts ");

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << infileList.Data() << endl;
        return 1;
    }

    // select categories here
    if (!loop.setInput("-*,+HParticleCandSim,+HGeantKine"))
    {
        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory *catParticle = loop.getCategory("HParticleCandSim");
    if (!catParticle)
    {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }
    HCategory *catGeant = loop.getCategory("HGeantKine");
    if (!catGeant)
    {
        std::cout << "No kineCat in input!" << std::endl;
        exit(1);
    }

    Int_t entries = loop.getEntries();
    if (nEvents < entries && nEvents >= 0)
        entries = nEvents;

    int primVertexBeforeDecayVertex = 0;
    int decayVertexBeforePrimVertex = 0;
    int primVertexInsideDecayVertex = 0;
    int decayVertexInsidePrimVertex = 0;

    // start of the event loop
    for (Int_t i = 1; i < nEvents; i++)
    {
        //----------break if last event is reached-------------
        if (loop.nextEvent(i) <= 0)
        {
            cout << " end recieved " << endl;
            break;
        } // last event reached
        HTool::printProgress(i, nEvents, 1, "Analysing evt# :");

        // for each event there are a number of tracks
        Int_t ntracks = catParticle->getEntries();

        std::vector<HRefitCand> protons, pions, kaons;
        std::vector<HParticleCandSim *> virtualCandLambdas, virtualCandProtons, virtualCandPions, virtualCandKaons;

        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim *cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack())
                continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed))
                continue;

            hGeantVertexXAll->Fill(cand->getGeantxVertex());
            hGeantVertexYAll->Fill(cand->getGeantyVertex());
            hGeantVertexZAll->Fill(cand->getGeantzVertex());

            HRefitCand candidate(cand);

            // select particles based on MC info
            // proton pdg==14, pion pdg==9, k+ pdg==11, lambda pdg==18
            // error values obtained from resoultion plots

            // The code below makes sure that all protons and pions in the analysis come from the Lambda decay
            // No combinatorial background

            /*if (cand->getGeantPID() == 14 && cand->getGeantParentPID() == 18) //Proton found
            {
                virtualCandProtons.push_back(cand);
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9 && cand->getGeantParentPID() == 18) // Pion found
            {
                virtualCandPions.push_back(cand);
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }*/

            if (cand->getGeantPID() == 14) //Proton found
            {
                virtualCandProtons.push_back(cand);
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9) // Pion found
            {
                virtualCandPions.push_back(cand);
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }
            else if (cand->getGeantPID() == 11) // Kaon found
            {
                virtualCandKaons.push_back(cand);
                double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                                   1.404, 2.723};
                FillData(cand, candidate, errors, 493.7);
                kaons.push_back(candidate);
            }
            if (cand->getGeantPID() == 18) //Lambda found
            {
                // There are probably no candidates for Lambda due to track requirement
                virtualCandLambdas.push_back(cand);
            }
            else
                continue;
        } // end track loop

        // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------

        //std::cout << "Event number: "  << i << std::endl;
        //std::cout << " " << std::endl;
        //std::cout << "Number of Kaons: " << kaons.size() << std::endl;
        //std::cout << "Number of protons: " << protons.size() << std::endl;
        //std::cout << "Number of Pions: " << pions.size() << std::endl;

        for (size_t s = 0; s < virtualCandLambdas.size(); s++)
        {

            HParticleCandSim *virtualCandLambda = virtualCandLambdas[s];

            hGeantVertexXLambda->Fill(virtualCandLambda->getGeantxVertex());
            hGeantVertexYLambda->Fill(virtualCandLambda->getGeantyVertex());
            hGeantVertexZLambda->Fill(virtualCandLambda->getGeantzVertex());
        }

        for (size_t n = 0; n < protons.size(); n++)
        {

            HRefitCand cand1 = protons[n];

            for (size_t m = 0; m < pions.size(); m++)
            {
                HRefitCand cand2 = pions[m];

                for (size_t p = 0; p < kaons.size(); p++)
                {
                    HRefitCand cand3 = kaons[p];

                    std::vector<HRefitCand> cands;
                    cands.clear();
                    cands.push_back(cand1);
                    cands.push_back(cand2);
                    cands.push_back(cand3);

                    // Initiate the vertex fitter
                    HVertexFinder vtxFinder(cands);

                    //vtxFitter.setVerbosity(0);
                    // Find the vertex

                    TVector3 vertex = vtxFinder.findVertex(cands);
                    TVector3 primaryVertex = vtxFinder.getPrimaryVertex();
                    //TVector3 primaryVertex = vtxFinder.findPrimaryVertex(cands);

                    hVertexXPreFit->Fill(vertex.X());
                    hVertexYPreFit->Fill(vertex.Y());
                    hVertexZPreFit->Fill(vertex.Z());

                    hPrimaryVertexXPreFit->Fill(primaryVertex.X());
                    hPrimaryVertexYPreFit->Fill(primaryVertex.Y());
                    hPrimaryVertexZPreFit->Fill(primaryVertex.Z());

                    double distBetweenVertices = vtxFinder.getDistBetweenVertices();
                    hDistPrimToDecayVertex->Fill(distBetweenVertices);

                    if (vtxFinder.isPrimVertexBeforeDecayVertex() == true)
                    {
                        primVertexBeforeDecayVertex++;
                    }
                    if (vtxFinder.isPrimVertexBeforeDecayVertex() == false)
                    {
                        decayVertexBeforePrimVertex++;
                    }
                    if (vtxFinder.isPrimVertexInsideDecayVertex() == false)
                    {
                        decayVertexInsidePrimVertex++;
                    }
                    if (vtxFinder.isPrimVertexInsideDecayVertex() == true)
                    {
                        primVertexInsideDecayVertex++;
                    }

                    // Create a vector of candidates to originate from primary vertex
                    std::vector<HRefitCand> primCands;
                    primCands.clear();
                    primCands.push_back(cand1);
                    primCands.push_back(cand3);

                    // Perform fitting of primary vertex
                    HVertexFitter vtxFitterPrimCands(primCands);
                    vtxFitterPrimCands.addVertexConstraint();
                    vtxFitterPrimCands.fit();
                    double probPrim;
                    probPrim = vtxFitterPrimCands.getProb();
                    std::cout << "Probability Prim vtx: " << probPrim << std::endl;

                    // Create a vector of candidates to originate from decay vertex
                    std::vector<HRefitCand> secCands;
                    secCands.clear();
                    secCands.push_back(cand1);
                    secCands.push_back(cand2);

                    // Perform fitting of secondary vertex
                    HVertexFitter vtxFitterSecCands(secCands);
                    vtxFitterSecCands.addVertexConstraint();
                    vtxFitterSecCands.fit();
                    double probSec;
                    probSec = vtxFitterSecCands.getProb();
                    std::cout << "Probability Sec vtx: " << probSec << std::endl;
                }
            }
        }

    } // end of the events loop

    std::cout << "Number of prim vertices before the decay vertex: " << primVertexBeforeDecayVertex << ", and after: " << decayVertexBeforePrimVertex << std::endl;
    std::cout << "Number of prim vertices inside the decay vertex: " << primVertexInsideDecayVertex << ", and outside: " << decayVertexInsidePrimVertex << std::endl;

    // write histograms to the output file
    outfile->cd();
    h01->Write();
    h02->Write();
    h03->Write();
    h04->Write();
    h05->Write();
    h06->Write();
    h07->Write();
    h08->Write();
    h09->Write();
    h10->Write();
    h11->Write();
    h12->Write();
    h13->Write();
    h14->Write();
    h15->Write();
    h16->Write();
    h17->Write();
    h18->Write();
    h19->Write();
    h20->Write();
    h21->Write();
    h22->Write();
    h23->Write();
    h24->Write();

    hVertexXPreFit->Write();
    hVertexYPreFit->Write();
    hVertexZPreFit->Write();
    hPrimaryVertexXPreFit->Write();
    hPrimaryVertexYPreFit->Write();
    hPrimaryVertexZPreFit->Write();

    hVertexPreFit->Write();
    hDistanceBetweenProtonAndPionPreFit->Write();
    hDistanceToVertexProtonPreFit->Write();
    hDistanceToVertexPionPreFit->Write();
    hDistanceToOriginProtonPreFit->Write();
    hDistanceToOriginPionPreFit->Write();
    hDistanceToOriginSumPreFit->Write();

    hVertexXPostFit->Write();
    hVertexYPostFit->Write();
    hVertexZPostFit->Write();
    hVertexPostFit->Write();

    hDistanceBetweenProtonAndPionPostFit->Write();
    hDistanceToVertexProtonPostFit->Write();
    hDistanceToVertexPionPostFit->Write();

    hGeantVertexXAll->Write();
    hGeantVertexYAll->Write();
    hGeantVertexZAll->Write();

    hGeantVertexXLambda->Write();
    hGeantVertexYLambda->Write();
    hGeantVertexZLambda->Write();

    hGeantVertexXCand1->Write();
    hGeantVertexYCand1->Write();
    hGeantVertexZCand1->Write();

    hGeantVertexXCand2->Write();
    hGeantVertexYCand2->Write();
    hGeantVertexZCand2->Write();

    hVertexXDiff->Write();
    hVertexYDiff->Write();
    hVertexZDiff->Write();

    hGeantTotMomentumProtons->Write();
    hGeantXMomentumProtons->Write();
    hGeantYMomentumProtons->Write();
    hGeantZMomentumProtons->Write();
    hGeantTotMomentumPions->Write();
    hGeantXMomentumPions->Write();
    hGeantYMomentumPions->Write();
    hGeantZMomentumPions->Write();

    hRecoRProtons->Write();
    hRecoZProtons->Write();
    hRecoThetaProtons->Write();
    hRecoPhiProtons->Write();
    hRecoMomentumProtons->Write();
    hRecoBetaProtons->Write();

    hRecoRPions->Write();
    hRecoZPions->Write();
    hRecoThetaPions->Write();
    hRecoPhiPions->Write();
    hRecoMomentumPions->Write();
    hRecoBetaPions->Write();

    h25->Write();
    h26->Write();
    h27->Write();
    h28->Write();
    h29->Write();
    h30->Write();
    h31->Write();
    h32->Write();
    h33->Write();

    hMomLambda->Write();
    hRecoThetaLambda->Write();
    hRecoPhiLambda->Write();

    hErrorRLambda->Write();
    hErrorZLambda->Write();
    hErrorThetaLambda->Write();
    hErrorPhiLambda->Write();

    hMomLambdaCut->Write();
    hRecoThetaLambdaCut->Write();
    hRecoPhiLambdaCut->Write();

    hErrorRLambdaCut->Write();
    hErrorZLambdaCut->Write();
    hErrorThetaLambdaCut->Write();
    hErrorPhiLambdaCut->Write();

    hDistPrimToDecayVertex->Write();

    outfile->Close();

    return 0;
} // end of the macro
