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

Int_t analysisVertexFitter(TString infileList = "pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents = 100000)
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

            if (cand->getGeantPID() == 14 && cand->getGeantParentPID() == 18) //Proton found
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
            }

            /*if (cand->getGeantPID() == 14) //Proton found
            {
                virtualCandProtons.push_back(cand);
                double errors[] = {1.469 * 1e-5*(1-0.05), 2.410 * 1e-3*(1-0.05), 5.895 * 1e-3*(1-0.05),
                                   1.188*(1-0.05), 2.652*(1-0.05)};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9) // Pion found
            {
                virtualCandPions.push_back(cand);
                double errors[] = {5.959 * 1e-5*(1-0.05), 9.316 * 1e-3*(1-0.05), 1.991 * 1e-2*(1-0.05),
                                   4.006*(1-0.05), 7.629*(1-0.05)};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }*/
            else if (cand->getGeantPID() == 11) // Kaon found
            {
                virtualCandKaons.push_back(cand);
                double errors[] = {1.947 * 1e-5*(1-0.05), 2.296 * 1e-3*(1-0.05), 6.312 * 1e-3*(1-0.05),
                                   1.404*(1-0.05), 2.723*(1-0.05)};
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
            HParticleCandSim *virtualCand1 = virtualCandProtons[n];

            for (size_t m = 0; m < pions.size(); m++)
            {
                HRefitCand cand2 = pions[m];
                HParticleCandSim *virtualCand2 = virtualCandPions[m];

                hRecoRProtons->Fill(cand1.getR());
                hRecoZProtons->Fill(cand1.getZ());
                hRecoThetaProtons->Fill(cand1.Theta());
                hRecoPhiProtons->Fill(cand1.Phi());
                hRecoMomentumProtons->Fill(cand1.P());
                //hRecoBetaProtons->Fill(cand1.getBeta());

                hRecoRPions->Fill(cand2.getR());
                hRecoZPions->Fill(cand2.getZ());
                hRecoThetaPions->Fill(cand2.Theta());
                hRecoPhiPions->Fill(cand2.Phi());
                hRecoMomentumPions->Fill(cand2.P());
                //hRecoBetaPions->Fill(cand2.getBeta());

                std::vector<HRefitCand> cands;
                cands.clear();
                cands.push_back(cand1);
                cands.push_back(cand2);

                // Initiate the vertex fitter
                HVertexFitter vtxFitter(cands);

                //vtxFitter.setVerbosity(0);
                // Find the vertex
                TVector3 vertex = vtxFitter.findVertex(cands);
                
		HVirtualCand lambdaCand;
		lambdaCand  = vtxFitter.getLambdaCandidate();

		std::cout << "--------------------" << std::endl;
		std::cout << "Lambda candidate: " << std::endl;
		std::cout << "Estimated Momentum: " << lambdaCand.getMomentum()  << std::endl;
		std::cout << "Theta: " << lambdaCand.getTheta() << std::endl;
		std::cout << "Phi: " << lambdaCand.getPhi() << std::endl;
                std::cout << "Theta: " << lambdaCand.Theta() << std::endl;
                std::cout << "Phi: " << lambdaCand.Phi() << std::endl;
                std::cout << "R: " << lambdaCand.getR() << std::endl;
                std::cout << "Z: " << lambdaCand.getZ() << std::endl;
                std::cout << "--------------------" << std::endl;
		
		TMatrixD lambdaCov(5, 5);
		lambdaCov=vtxFitter.getCovarianceMatrixLambda();
		lambdaCov.Print();
	    	
		double deg2rad = TMath::DegToRad();
                double rad2deg = TMath::RadToDeg();

		//HVirtualCand* lambdaCandPtr;
		//lambdaCandPtr = Â&lambdaCand;

		HRefitCand lambdaCandRefit(&lambdaCand);
		lambdaCandRefit.SetXYZM(lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                        std::cos(lambdaCand.getPhi() * deg2rad),
                    lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                        std::sin(lambdaCand.getPhi() * deg2rad),
                    lambdaCand.getMomentum() * std::cos(lambdaCand.getTheta() * deg2rad),
                    1115.683);
    		lambdaCandRefit.setR(lambdaCand.getR());
    		lambdaCandRefit.setZ(lambdaCand.getZ());
    		lambdaCandRefit.setCovariance(lambdaCov);

		std::cout << "--------------------" << std::endl;
                std::cout << "Lambda candidate Refit: " << std::endl;
                std::cout << "Estimated Momentum: " << lambdaCandRefit.P()  << std::endl;
                std::cout << "Theta: " << lambdaCandRefit.Theta() << " , in degrees: " << lambdaCandRefit.Theta()*rad2deg << std::endl;
                std::cout << "Phi: " << lambdaCandRefit.Phi() << " , in degrees: " << lambdaCandRefit.Phi()*rad2deg << std::endl;
                std::cout << "R: " << lambdaCandRefit.getR() << std::endl;
                std::cout << "Z: " << lambdaCandRefit.getZ() << std::endl;
                std::cout << "--------------------" << std::endl;

		TMatrixD lambdaRefitCov(5,5);
                lambdaRefitCov=lambdaCandRefit.getCovariance();
		lambdaRefitCov.Print();
		
		double distProtonPion = vtxFitter.getDistanceBetweenFittedParticles();

                hDistanceToVertexProtonPreFit->Fill(vtxFitter.getDistanceFirstParticleVertex());
                hDistanceToVertexPionPreFit->Fill(vtxFitter.getDistanceSecondParticleVertex());

                hDistanceToOriginProtonPreFit->Fill(vtxFitter.getDistanceFirstParticleOrigin());
                hDistanceToOriginPionPreFit->Fill(vtxFitter.getDistanceSecondParticleOrigin());

                hDistanceToOriginSumPreFit->Fill(vtxFitter.getDistanceFirstParticleOrigin() + vtxFitter.getDistanceSecondParticleOrigin());

                // Write out the components of the vertex vector
                //std::cout << "Vertex pos: x=" << vertex.X() << ", y=" << vertex.Y() << ", z=" << vertex.Z() << std::endl;

                hVertexXPreFit->Fill(vertex.X());
                hVertexYPreFit->Fill(vertex.Y());
                hVertexZPreFit->Fill(vertex.Z());

                hGeantVertexXCand1->Fill(virtualCand1->getGeantxVertex());
                hGeantVertexYCand1->Fill(virtualCand1->getGeantyVertex());
                hGeantVertexZCand1->Fill(virtualCand1->getGeantzVertex());

                hGeantVertexXCand2->Fill(virtualCand2->getGeantxVertex());
                hGeantVertexYCand2->Fill(virtualCand2->getGeantyVertex());
                hGeantVertexZCand2->Fill(virtualCand2->getGeantzVertex());

                hVertexXDiff->Fill(vertex.X() - virtualCand1->getGeantxVertex());
                hVertexYDiff->Fill(vertex.Y() - virtualCand1->getGeantyVertex());
                hVertexZDiff->Fill(vertex.Z() - virtualCand1->getGeantzVertex());

                hGeantTotMomentumProtons->Fill(virtualCand1->getGeantTotalMom());
                hGeantXMomentumProtons->Fill(virtualCand1->getGeantxMom());
                hGeantYMomentumProtons->Fill(virtualCand1->getGeantyMom());
                hGeantZMomentumProtons->Fill(virtualCand1->getGeantzMom());
                hGeantTotMomentumPions->Fill(virtualCand2->getGeantTotalMom());
                hGeantXMomentumPions->Fill(virtualCand2->getGeantxMom());
                hGeantYMomentumPions->Fill(virtualCand2->getGeantyMom());
                hGeantZMomentumPions->Fill(virtualCand2->getGeantzMom());

                double R = sqrt((vertex.X() * vertex.X()) + (vertex.Y() * vertex.Y()));
                hVertexPreFit->Fill(vertex.Z(), R);

                hDistanceBetweenProtonAndPionPreFit->Fill(distProtonPion);
                //hDistanceToVertexProtonPreFit->Fill();
                //hDistanceToVertexPionPreFit->Fill();

                std::vector<HRefitCand> newCands = vtxFitter.UpdateTrackParameters(cands, vertex);

                HRefitCand newCand1=newCands[0];
                newCand1.setZ(vertex.Z());
                if(std::cos(newCand1.Phi()+ TMath::PiOver2())!=0){
                double newR1 = vertex.X()/std::cos(newCand1.Phi()+ TMath::PiOver2());
                newCand1.setR(newR1);
                }

                HRefitCand newCand2=newCands[1];
                newCand2.setZ(vertex.Z());
                if(std::cos(newCand1.Phi()+ TMath::PiOver2())!=0){
                double newR2 = vertex.X()/std::cos(newCand2.Phi()+ TMath::PiOver2());
                newCand2.setR(newR2);
                }

                std::vector<HRefitCand> updatedCands;
                updatedCands.clear();
                updatedCands.push_back(newCand1);
                updatedCands.push_back(newCand2);

                // Fit to secondary vertex
                //HVertexFitter vtxFitterNew(newCands);
                // Fit to primary vertex
                HVertexFitter vtxFitterNew(updatedCands);

                vtxFitterNew.setPhiOriginal(cand1.Phi(), cand2.Phi());
                vtxFitterNew.setLearningRate(1);
                vtxFitterNew.setNumberOfIterations(5);

                vtxFitterNew.fit();

                if (vtxFitterNew.isConverged() == true)
                {
                    h10->Fill(1);
                }
                else if (vtxFitterNew.isConverged() == false)
                {
                    h10->Fill(0);

                    // Filling the Pull histos for proton
                    h25->Fill(vtxFitterNew.getPull(0));
                    h26->Fill(vtxFitterNew.getPull(1));
                    h27->Fill(vtxFitterNew.getPull(2));
                    h28->Fill(vtxFitterNew.getPull(3));
                    h29->Fill(vtxFitterNew.getPull(4));
                    // Filling the Pull histos for pion
                    h30->Fill(vtxFitterNew.getPull(6));
                    h31->Fill(vtxFitterNew.getPull(7));
                    h32->Fill(vtxFitterNew.getPull(8));
                    h33->Fill(vtxFitterNew.getPull(9));


                }

                h11->Fill(vtxFitterNew.getIteration());

                h02->Fill(vtxFitterNew.getChi2());
                h03->Fill(vtxFitterNew.getProb());

                //get Pull example (1/P for the fitted proton)
                h05->Fill(vtxFitterNew.getPull(0));

                //get Pull example (theta for the fitted proton)
                h06->Fill(vtxFitterNew.getPull(1));

                //get Pull example (phi for the fitted proton)
                h07->Fill(vtxFitterNew.getPull(2));

                //get Pull example (R for the fitted proton)
                h08->Fill(vtxFitterNew.getPull(3));

                //get Pull example (Z for the fitted proton)
                h09->Fill(vtxFitterNew.getPull(4));

                h17->Fill(vtxFitterNew.getPull(6));

                h18->Fill(vtxFitterNew.getPull(7));

                h19->Fill(vtxFitterNew.getPull(8));

                h20->Fill(vtxFitterNew.getPull(9));

                if (vtxFitterNew.getProb() > 0.01)
                {

                    // get Pull example (1/P for the fitted proton)
                    h12->Fill(vtxFitterNew.getPull(0));

                    // get Pull example (theta for the fitted proton)
                    h13->Fill(vtxFitterNew.getPull(1));

                    // get Pull example (phi for the fitted proton)
                    h14->Fill(vtxFitterNew.getPull(2));

                    // get Pull example (R for the fitted proton)
                    h15->Fill(vtxFitterNew.getPull(3));

                    // get Pull example (Z for the fitted proton)
                    h16->Fill(vtxFitterNew.getPull(4));

                    h21->Fill(vtxFitterNew.getPull(6));

                    h22->Fill(vtxFitterNew.getPull(7));

                    h23->Fill(vtxFitterNew.getPull(8));

                    h24->Fill(vtxFitterNew.getPull(9));
                }
            }
        }

        /*       // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------
        for (size_t n = 0; n < protons.size(); n++)
        {
            HRefitCand cand1 = protons[n];
            for (size_t m = 0; m < pions.size(); m++)
            {
                HRefitCand cand2 = pions[m];
                // mass prefit
                h01->Fill((cand1 + cand2).M());

                // ---------------------------------------------------------------------------------
                // begin kinfit here
                // ---------------------------------------------------------------------------------
                std::vector<HRefitCand> cands;
                cands.push_back(cand1);
                cands.push_back(cand2);

                HFitter fitter(cands);
                fitter.addMassConstraint(1115.68);
                double learningRate=1.0;
		int numIter=1;
		fitter.fit(learningRate, numIter);

                // get fitted objects fittedcand1 and fittedcand2
                HRefitCand fcand1 = fitter.getDaughter(0); // proton
                HRefitCand fcand2 = fitter.getDaughter(1); // pion

                h02->Fill(fitter.getChi2());
                h03->Fill(fitter.getProb());
                h04->Fill((fcand1 + fcand2).M());

                 //get Pull example (1/P for the fitted proton)
                h05->Fill(fitter.getPull(0));
		
		//get Pull example (theta for the fitted proton)
		h06->Fill(fitter.getPull(1));
		
		//get Pull example (phi for the fitted proton)
                h07->Fill(fitter.getPull(2));
		
		//get Pull example (R for the fitted proton)
                h08->Fill(fitter.getPull(3));
		
		//get Pull example (Z for the fitted proton)
                h09->Fill(fitter.getPull(4));
		
		if(fitter.isConverged()==true){
		h10->Fill((fcand1 + fcand2).M());
		}
		
		if(fitter.getProb()>0.01){
		               
		h17->Fill((fcand1 + fcand2).M());
		
		// get Pull example (1/P for the fitted proton)
                h12->Fill(fitter.getPull(0));
		
		// get Pull example (theta for the fitted proton)
		h13->Fill(fitter.getPull(1));
		
		// get Pull example (phi for the fitted proton)
                h14->Fill(fitter.getPull(2));
		
		// get Pull example (R for the fitted proton)
                h15->Fill(fitter.getPull(3));
		
		// get Pull example (Z for the fitted proton)
                h16->Fill(fitter.getPull(4));
		
		}
		h11->Fill(fitter.getIteration());
		
                // ---------------------------------------------------------------------------------
            }
        } */
        // -----------------------------------------------------------------------
        // -----------------------------------------------------------------------
        // Apply the kinfit with vertex constraint to proton and kaon
        // -----------------------------------------------------------------------
        /*        for (size_t n = 0; n < protons.size(); n++)
        {
            HRefitCand cand1 = protons[n];
            for (size_t m = 0; m < kaons.size(); m++)
            {
                HRefitCand cand2 = kaons[m];
                // mass prefit
                h01->Fill((cand1 + cand2).M());

                // ---------------------------------------------------------------------------------
                // begin kinfit here
                // ---------------------------------------------------------------------------------
                std::vector<HRefitCand> cands;
                cands.push_back(cand1);
                cands.push_back(cand2);

                HVertexFitter fitter(cands);
		fitter.addVtxConstraint();
                double learningRate=1;
		int numIter=1;
		
		// Using my modifications to KinFit
		//fitter.fit(learningRate, numIter);
		fitter.fit();
                
		// get fitted objects fittedcand1 and fittedcand2
                HRefitCand fcand1 = fitter.getDaughter(0); // proton
                HRefitCand fcand2 = fitter.getDaughter(1); // kaon

                h02->Fill(fitter.getChi2());
                h03->Fill(fitter.getProb());
                h04->Fill((fcand1 + fcand2).M());

                // get Pull example (1/P for the fitted proton)
                h05->Fill(fitter.getPull(0));
		
		// get Pull example (theta for the fitted proton)
		h06->Fill(fitter.getPull(1));
		
		// get Pull example (phi for the fitted proton)
                h07->Fill(fitter.getPull(2));
		
		// get Pull example (R for the fitted proton)
                h08->Fill(fitter.getPull(3));
		
		// get Pull example (Z for the fitted proton)
                h09->Fill(fitter.getPull(4));
		
		if(fitter.getProb()>0.01){
		
		// get Pull example (1/P for the fitted proton)
                h12->Fill(fitter.getPull(0));
		
		// get Pull example (theta for the fitted proton)
		h13->Fill(fitter.getPull(1));
		
		// get Pull example (phi for the fitted proton)
                h14->Fill(fitter.getPull(2));
		
		// get Pull example (R for the fitted proton)
                h15->Fill(fitter.getPull(3));
		
		// get Pull example (Z for the fitted proton)
                h16->Fill(fitter.getPull(4));
		
		}
		
		if(fitter.isConverged()==true){
		h10->Fill((fcand1 + fcand2).M());
		}
		
		h11->Fill(fitter.getIteration());
		
                // ---------------------------------------------------------------------------------
            }
        } */

    } // end of the events loop

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

    outfile->Close();

    return 0;
} // end of the macro
