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

Int_t analysis(TString infileList = "/lustre/hades/user/jregina/DstTest/pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents = -1)
{

    // My input: pp_pKlambda_100000evts1_dst_apr12.root

    TStopwatch timer;
    timer.Start();

    //Momentum dependent uncertainty estimation input
    TFile *momErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepMomErrors_allParticles_errFunc.root", "read");
    TF1 *momErrP = (TF1 *)momErr->Get("f_pP");
    TF1 *momErrPi = (TF1 *)momErr->Get("f_pPi");
    TF1 *momErrK = (TF1 *)momErr->Get("f_pK");

    TFile *thtErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_allParticles_errFunc.root", "read");
    TFile *thtErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_PK_errFunc.root", "read");
    TF1 *thtErrP = (TF1 *)thtErr->Get("f_thtP");
    TF1 *thtErrPi = (TF1 *)thtErr->Get("f_thtPi");
    TF1 *thtErrK = (TF1 *)thtErr_PK->Get("f_thtK");

    TFile *phiErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PPi_errFunc.root", "read"); //if 4500 in name this is just a name issue
    TFile *phiErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PK_errFunc.root", "read");
    TF1 *phiErrP = (TF1 *)phiErr->Get("f_phiP");
    TF1 *phiErrPi = (TF1 *)phiErr->Get("f_phiPi");
    TF1 *phiErrK = (TF1 *)phiErr_PK->Get("f_phiK");

    TFile *RZErr_PPi = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PPi_errFunc.root", "read");
    TFile *RZErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PK_errFunc.root", "read");
    TF1 *RErrP = (TF1 *)RZErr_PPi->Get("f_RP");
    TF1 *RErrPi = (TF1 *)RZErr_PPi->Get("f_RPi");
    TF1 *RErrK = (TF1 *)RZErr_PK->Get("f_RK");
    TF1 *ZErrP = (TF1 *)RZErr_PPi->Get("f_ZP");
    TF1 *ZErrPi = (TF1 *)RZErr_PPi->Get("f_ZPi");
    TF1 *ZErrK = (TF1 *)RZErr_PK->Get("f_ZK");

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile *outfile = new TFile("testvertexfit.root", "recreate");

    TH1F *h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F *h02 = new TH1F("hChi2", "", 100, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" Counts ");
    TH1F *h021 = (TH1F *)h02->Clone("hChi2_3c");
    TH1F *h023 = (TH1F *)h02->Clone("hChi2_3cConv");

    TH1F *h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" Counts ");
    TH1F *h031 = (TH1F *)h03->Clone("hPChi2_3c");
    TH1F *h033 = (TH1F *)h03->Clone("hPChi2_3cConv");

    TH1F *h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" Events ");
    TH1F *h040 = (TH1F *)h04->Clone("hLambdaMassPostFitCut");
    TH1F *h041 = (TH1F *)h04->Clone("hLambdaMassPostFit_3c");
    TH1F *h044 = (TH1F *)h04->Clone("hLambdaMassPostFit_3cConvCut");

    TH1F *h05 = new TH1F("hPullPInv", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" Counts ");
    TH1F *h050 = (TH1F *)h05->Clone("hPullPInvCut");
    TH1F *h051 = (TH1F *)h05->Clone("hPullPInv_3c");
    TH1F *h054 = (TH1F *)h05->Clone("hPullPInv_3cConvCut");
    /*
    TH1F* h06 = new TH1F("hTotMomPreFit", "", 100, 3800, 4800);
    h06->SetXTitle(" p [MeV/c]");
    h06->SetYTitle(" events ");

    TH1F* h07 = new TH1F("hTotMomPostFitters", "", 100, 3800, 4800);
    h07->SetXTitle(" p [MeV/c]");
    h07->SetYTitle(" events ");
    TH1F *h074 = (TH1F*)h07->Clone("hTotMomPostFit_3cConvCut");
*/
    TH1F *h10 = new TH1F("hSuccessfulConvergence", "", 2, 0, 2);
    h10->SetXTitle("Fit has converged");
    h10->SetYTitle(" Counts ");

    TH1F *h08 = new TH1F("hIterations", "", 10, 0, 10);
    h08->SetXTitle("Number of Iterations");
    h08->SetYTitle(" Counts ");
    TH1F *h081 = (TH1F *)h08->Clone("hNIterations_3c");
    TH1F *h084 = (TH1F *)h08->Clone("hNIterations_3cCut");

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

    // Vertex prim and decay info
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
    if (nEvents > entries || nEvents < 0)
        nEvents = entries;

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

            if (cand->getGeantPID() == 14) //Proton found
            {
                virtualCandProtons.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                //                  1.188, 2.652};   //rough error estimates
                double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                                   RErrP->Eval(mom), ZErrP->Eval(mom)}; //momentum dependent error estimates
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9) // Pion found
            {
                virtualCandPions.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                //                  4.006, 7.629};   //rough error estimates
                double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                                   RErrPi->Eval(mom), ZErrPi->Eval(mom)}; //momentum dependent error estimates
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }
            else if (cand->getGeantPID() == 11) // Kaon found
            {
                virtualCandKaons.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                //                 1.404, 2.723};   //rough error estimates
                double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                                   RErrK->Eval(mom), ZErrK->Eval(mom)}; //momentum dependent error estimates
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

        primVtxFound = false;
        secVtxFound = false;

        //if (protons.size() < 2)
        //  continue;

        TVector3 primaryVertex, decayVertex;

        double probPrim = -99999, probSec = -99999;

        bool bestDecayVertexFound = false;
        bool bestPrimVertexFound = false;

        TVector3 decayVertex_Temp, primVertex_Temp;
        TVector3 decayVertex_TempUpdated, tempVertex_TempUpdated;
        TVector3 decayVertexBestFit, primVertexBestFit;

        double probDecayVertex_Temp = -1;
        double probPrimVertex_Temp = -1;
        double probDecayVertex_TempUpdated = -1;
        double probPrimVertex_TempUpdated = -1;

        double bestDiffXDecay = -1;
        double bestDiffYDecay = -1;
        double bestDiffZDecay = -1;
        double bestDiffXPrim = -1;
        double bestDiffYPrim = -1;
        double bestDiffZPrim = -1;

        double indexPrimaryProton;
        double indexDecayProton;

        std::vector<HRefitCand> cands3c;
        cands3c.clear();

        // Brief description
        // Looping over all found protons in one event.
        // For each found proton, all found pions and kaons are looped over respectively
        // A vertex is found and fitted for all combinations of protons and pions
        // A vertex found from a proton-kaon combination is labelled as a primary vertex
        // A vertex found from a proton-pion combination is labelled as a decay vertex
        // If several primary vertices are found in one event, the primary vertex with the best (highest) fit probability is kept
        // If several decay vertices are found in one event, the decay vertex with the best (highest) fit probability is kept
        // If both a primary vertex and a decay vertex were found in the same event, a neutral mother candidate is created with 
        // parameters calculated from the two vertices
        // A 3C fit is performed to the decay products and the neutral mother candidate utilizing momentum conservation

        // Proton loop
        for (size_t n = 0; n < protons.size(); n++)
        {

            // Proton
            HRefitCand cand1 = protons[n];

            //HParticleCandSim *virtualCand1 = virtualCandProtons[n];

            // Pion loop
            for (size_t m = 0; m < pions.size(); m++)
            {
                
                // Pion
                HRefitCand cand2 = pions[m];
                
                // Object for obtaining the Geant information and 
                // calculating the resolution
                HParticleCandSim *virtualCand2 = virtualCandPions[m];

                // A decay vertex has been found
                secVtxFound = true;
                
                // Vector of HRefitCand to pass to the vertex fit
                std::vector<HRefitCand> candsSec;
                candsSec.clear();
                candsSec.push_back(cand1);
                candsSec.push_back(cand2);

                // Calculating the invariant mass before the fit
                TLorentzVector lambda = cand1 + cand2;
                h045->Fill(lambda.M());

                // Initiating the vertex finder
                HVertexFinder *vtxFinderSec = new HVertexFinder();
                
                // Finding the decay vertex
                decayVertex = vtxFinderSec->findVertex(candsSec);

                // Filling histograms with vertex positions and distances between particles and particles-decay vertex
                hVertexXPreFit->Fill(decayVertex.X());
                hVertexYPreFit->Fill(decayVertex.Y());
                hVertexZPreFit->Fill(decayVertex.Z());
                hDistanceToVertexProtonPreFit->Fill(vtxFinderSec->getDistanceFirstParticleVertex());
                hDistanceToVertexPionPreFit->Fill(vtxFinderSec->getDistanceSecondParticleVertex());
                hDistanceBetweenProtonAndPionPreFit->Fill(vtxFinderSec->getDistanceBetweenFittedParticles());

                // Perform fitting of secondary vertex
                HVertexFitter vtxFitterSecCands(candsSec);
                
                // Adding the vertex constraint
                vtxFitterSecCands.addVertexConstraint();
                
                // Performing the fit
                vtxFitterSecCands.fit();
                
                // Filling histograms with the chi2 and probability
                h02SecondaryVtx->Fill(vtxFitterSecCands.getChi2());
                h03SecondaryVtx->Fill(vtxFitterSecCands.getProb());

                // Example Pull for 1/p of proton
                h06SecondaryVtx->Fill(vtxFitterSecCands.getPull(1));

                probSec = vtxFitterSecCands.getProb();

                // If several decay vertices were found, the one with the largest probability is kept
                if (probSec > probDecayVertex_Temp)
                {

                    // Updating variables that are used if a best decay vertex was found
                    bestDecayVertexFound = true;
                    probDecayVertex_Temp = probSec;
                    decayVertexBestFit = decayVertex;
                    bestDiffXDecay = decayVertex.X() - virtualCand2->getGeantxVertex();
                    bestDiffYDecay = decayVertex.Y() - virtualCand2->getGeantyVertex();
                    bestDiffZDecay = decayVertex.Z() - virtualCand2->getGeantzVertex();

                    // Index of proton to be compared with the index of the proton that was used to find the primary vertex
                    indexDecayProton = n;

                    // Vector with charged particles to pass to the 3C fit
                    cands3c.clear();

                    // Get the proton daughter and pass it to the 3C fit later
                    cands3c.push_back(vtxFitterSecCands.getDaughter(0));

                    // Get the pion daughter and pass it to the 3C fit later
                    cands3c.push_back(vtxFitterSecCands.getDaughter(1));
                }

                // Fill histogram with the number of iterations before the fit converged
                h11SecVtx->Fill(vtxFitterSecCands.getIteration());

                if (probSec > 0)
                {
                    // Example Pull for 1/p of proton
                    h13SecondaryVtx->Fill(vtxFitterSecCands.getPull(1));

                    // Vertex reolution plots
                    hVertexXDiff_ProbCut->Fill(decayVertex.X() - virtualCand2->getGeantxVertex());
                    hVertexYDiff_ProbCut->Fill(decayVertex.Y() - virtualCand2->getGeantyVertex());
                    hVertexZDiff_ProbCut->Fill(decayVertex.Z() - virtualCand2->getGeantzVertex());

                    // Create a neutral mother candidate from the fitted objects
                    TLorentzVector lambdaVtxFit = vtxFitterSecCands.getDaughter(0) + vtxFitterSecCands.getDaughter(1);

                    // Checking how the mass histogram looks after the decay vertex fit
                    h040->Fill(lambdaVtxFit.M());
                }

                // Vertex reolution plots
                hVertexXDiff->Fill(decayVertex.X() - virtualCand2->getGeantxVertex());
                hVertexYDiff->Fill(decayVertex.Y() - virtualCand2->getGeantyVertex());
                hVertexZDiff->Fill(decayVertex.Z() - virtualCand2->getGeantzVertex());

                // Reconstructed track parameters for the proton
                hRecoRProtons->Fill(cand1.getR());
                hRecoZProtons->Fill(cand1.getZ());
                hRecoThetaProtons->Fill(cand1.Theta());
                hRecoPhiProtons->Fill(cand1.Phi());
                hRecoMomentumProtons->Fill(cand1.P());

                // Reconstructed track parameters pions
                hRecoRPions->Fill(cand2.getR());
                hRecoZPions->Fill(cand2.getZ());
                hRecoThetaPions->Fill(cand2.Theta());
                hRecoPhiPions->Fill(cand2.Phi());
                hRecoMomentumPions->Fill(cand2.P());
            }

            // Kaon loop
            for (size_t p = 0; p < kaons.size(); p++)
            {

                // Kaons
                HRefitCand cand3 = kaons[p];

                // Object used for obtaining the Geant information for the resolution plots
                HParticleCandSim *virtualCand3 = virtualCandKaons[p];

                // A primary vertex was found in the event
                primVtxFound = true;

                // Vector with primary candidates
                std::vector<HRefitCand> candsPrim;
                candsPrim.clear();
                candsPrim.push_back(cand1);
                candsPrim.push_back(cand3);

                // Initiate the vertex finding
                HVertexFinder *vtxFinderPrim = new HVertexFinder();
                
                // Find a primary vertex from the current candidates
                primaryVertex = vtxFinderPrim->findVertex(candsPrim);

                // Vertex histograms
                hPrimaryVertexXPreFit->Fill(primaryVertex.X());
                hPrimaryVertexYPreFit->Fill(primaryVertex.Y());
                hPrimaryVertexZPreFit->Fill(primaryVertex.Z());

                // Histograms with distances between the particles and the particles-vertex
                hDistanceBetweenPrimaryProtonAndKaonPreFit->Fill(vtxFinderPrim->getDistanceBetweenFittedParticles());
                hDistanceToVertexPrimaryProtonPreFit->Fill(vtxFinderPrim->getDistanceFirstParticleVertex());
                hDistanceToVertexKaonPreFit->Fill(vtxFinderPrim->getDistanceSecondParticleVertex());

                // Reconstructed track parameters protons
                hRecoRPrimProtons->Fill(cand1.getR());
                hRecoZPrimProtons->Fill(cand1.getZ());
                hRecoThetaPrimProtons->Fill(cand1.Theta());
                hRecoPhiPrimProtons->Fill(cand1.Phi());
                hRecoMomentumPrimProtons->Fill(cand1.P());

                // Reconstructed track parameters kaons
                hRecoRKaons->Fill(cand3.getR());
                hRecoZKaons->Fill(cand3.getZ());
                hRecoThetaKaons->Fill(cand3.Theta());
                hRecoPhiKaons->Fill(cand3.Phi());
                hRecoMomentumKaons->Fill(cand3.P());

                // Initiating the vertex fitter
                HVertexFitter vtxFitterPrimCands(candsPrim);
                
                // Adding the vertex constraint
                vtxFitterPrimCands.addVertexConstraint();
                
                // Performing the fit
                vtxFitterPrimCands.fit();

                // FIlling histograms with the chi2 and probability
                h02->Fill(vtxFitterPrimCands.getChi2());
                h03->Fill(vtxFitterPrimCands.getProb());

                probPrim = vtxFitterPrimCands.getProb();

                if (probPrim > probPrimVertex_Temp)
                {

                    // Updating some variables that are saved for the primary vertex with the largest probability 
                    bestPrimVertexFound = true;
                    probPrimVertex_Temp = probPrim;
                    primVertexBestFit = primaryVertex;                    
                    bestDiffXPrim = primaryVertex.X() - virtualCand3->getGeantxVertex();
                    bestDiffYPrim = primaryVertex.Y() - virtualCand3->getGeantyVertex();
                    bestDiffZPrim = primaryVertex.Z() - virtualCand3->getGeantzVertex();
                    
                    // Index of the proton. To be compared later to the index of the proton from which the decay vertex was found
                    indexPrimaryProton = n;

                }

                // Fill histogram with the number of iterations before convergence
                h11->Fill(vtxFitterPrimCands.getIteration());

                // Example Pull for 1/p of proton
                h06->Fill(vtxFitterPrimCands.getPull(1));

                // Vertex resolutions
                hVertexXDiffPrim->Fill(primaryVertex.X() - virtualCand3->getGeantxVertex());
                hVertexYDiffPrim->Fill(primaryVertex.Y() - virtualCand3->getGeantyVertex());
                hVertexZDiffPrim->Fill(primaryVertex.Z() - virtualCand3->getGeantzVertex());

                if (probPrim > 0)
                {

                    // Vertex resolutions
                    hVertexXDiffPrim_ProbCut->Fill(primaryVertex.X() - virtualCand3->getGeantxVertex());
                    hVertexYDiffPrim_ProbCut->Fill(primaryVertex.Y() - virtualCand3->getGeantyVertex());
                    hVertexZDiffPrim_ProbCut->Fill(primaryVertex.Z() - virtualCand3->getGeantzVertex());

                    // Example Pull for 1/p of proton
                    h13->Fill(vtxFitterPrimCands.getPull(1));
                }
            }
        }

        // Testing if both aprimary and decay vertex wa found
        if (bestPrimVertexFound == true && bestDecayVertexFound == true)
        {
            // Making sure that two different protons were used to reconstruct the vertices
            if (indexDecayProton != indexPrimaryProton)
            {

                // Calculating the distance between the vertices
                double distPrimToDecayVertex = sqrt((decayVertexBestFit.X() - primVertexBestFit.X()) * (decayVertexBestFit.X() - primVertexBestFit.X()) + (decayVertexBestFit.Y() - primVertexBestFit.Y()) * (decayVertexBestFit.Y() - primVertexBestFit.Y()) + (decayVertexBestFit.Z() - primVertexBestFit.Z()) * (decayVertexBestFit.Z() - primVertexBestFit.Z()));

                hDistPrimToDecayVertex->Fill(distPrimToDecayVertex);

                // Vertex resolutions
                hVertexXDiff_BothVerticesFound->Fill(bestDiffXDecay);
                hVertexYDiff_BothVerticesFound->Fill(bestDiffYDecay);
                hVertexZDiff_BothVerticesFound->Fill(bestDiffZDecay);
                hVertexXDiffPrim_BothVerticesFound->Fill(bestDiffXPrim);
                hVertexYDiffPrim_BothVerticesFound->Fill(bestDiffYPrim);
                hVertexZDiffPrim_BothVerticesFound->Fill(bestDiffZPrim);

                double R_primaryVertex, R_decayVertex;

                R_primaryVertex = sqrt(primVertexBestFit.X() * primVertexBestFit.X() + primVertexBestFit.Y() * primVertexBestFit.Y());
                R_decayVertex = sqrt(decayVertexBestFit.X() * decayVertexBestFit.X() + decayVertexBestFit.Y() * decayVertexBestFit.Y());

                // Vertex histograms when a primary and a decay vertex was found in the same event
                hVertexXPostFit->Fill(decayVertexBestFit.X());
                hVertexYPostFit->Fill(decayVertexBestFit.Y());
                hVertexZPostFit->Fill(decayVertexBestFit.Z());
                hPrimVertexXPostFit->Fill(primVertexBestFit.X());
                hPrimVertexYPostFit->Fill(primVertexBestFit.Y());
                hPrimVertexZPostFit->Fill(primVertexBestFit.Z());

                // FInding the neutral mother candidate
                HNeutralCandFinder lambdaCandFinder(cands3c);
                
                // Make sure the primary vertex is used in the calculation of the neutral mother candidate
                lambdaCandFinder.setUsePrimaryVertexInNeutralMotherCalculation(true);
                
                // Calculated the properties of the neutral mother candidate
                lambdaCandFinder.setNeutralMotherCandFromPrimaryVtxInfo(primaryVertex, decayVertex);

                // Get the neutral mother candidate as a HVritualCand
                HVirtualCand lambdaCand = lambdaCandFinder.getNeutralMotherCandidate();

                // Creating an HRefitCand object that can be used in the 3C fit
                HRefitCand lambdaCandRefit(&lambdaCand);
                lambdaCandRefit.SetXYZM(lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                            std::cos(lambdaCand.getPhi() * deg2rad),
                                        lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                            std::sin(lambdaCand.getPhi() * deg2rad),
                                        lambdaCand.getMomentum() * std::cos(lambdaCand.getTheta() * deg2rad),
                                        1115.683);
                
                hMomLambda->Fill(lambdaCandRefit.P());
                
                // Setting the estimated covariance matrix for the neutral mother candidate
                TMatrixD lambdaCov(5, 5);
                lambdaCov = lambdaCandFinder.getCovarianceMatrixNeutralMother();
                lambdaCandRefit.setCovariance(lambdaCov);
            
                // Setting the parameters for the neutral mother candidate
                lambdaCandRefit.setR(lambdaCand.getR());
                lambdaCandRefit.setZ(lambdaCand.getZ());
                lambdaCandRefit.SetTheta(lambdaCand.getTheta() * deg2rad);
                lambdaCandRefit.SetPhi(lambdaCand.getPhi() * deg2rad);

                // Fill histograms with the neutral mother candidate parameters
                hRecoRLambda->Fill(lambdaCand.getR());
                hRecoZLambda->Fill(lambdaCand.getZ());
                hRecoThetaLambda->Fill(lambdaCand.getTheta() * deg2rad);
                hRecoPhiLambda->Fill(lambdaCand.getPhi() * deg2rad);

                // Take the lambda mass before the fit
                TLorentzVector lambdaCandBefore3C = cands3c[0] + cands3c[1];
                hLambdaMassBeforeFit->Fill(lambdaCandBefore3C.M());

                if (cands3c.size() == 2)
                {
                    // Adding the candidates to the 3C fit and the neutral mother candidate
                    HVertexFitter Fitter3c(cands3c, lambdaCandRefit);

                    // Setting the constraint and the maximum number of iterations
                    Fitter3c.add3Constraint();
                    Fitter3c.setNumberOfIterations(20);

                    // Performing the actual fit
                    Fitter3c.fit();
                    
                    
                    if (Fitter3c.isConverged())
                    {
                        
                        // Retrieving the mother candidate from after the fit
                        HRefitCand flambda = Fitter3c.getMother();

                        // Retrieving the daughter particles after the fit
                        HRefitCand cand13C = Fitter3c.getDaughter(0); // proton
                        HRefitCand cand23C = Fitter3c.getDaughter(1); // pion
                        
                        // Calculating the invariant mass of the daughter particles
                        TLorentzVector lambdaCand3C = cand13C + cand23C;
                        h04->Fill(lambdaCand3C.M());
                    }

                    // Fill histograms with the chi2 and probability of the 3C fit
                    h021->Fill(Fitter3c.getChi2());
                    h031->Fill(Fitter3c.getProb());
                }

                if (R_primaryVertex < R_decayVertex)
                {

                    primVertexInsideDecayVertex++;
                }
                else
                {
                    decayVertexInsidePrimVertex++;
                }

                // Fill histograms with the neutral mother candidate parameters if the decay vertex is located further downstream than the primary vertex
                if (primVertexBestFit.Z() < decayVertexBestFit.Z())
                {

                    hRecoRLambdaZCut->Fill(lambdaCand.getR());
                    hRecoZLambdaZCut->Fill(lambdaCand.getZ());
                    hRecoThetaLambdaZCut->Fill(lambdaCand.getTheta() * deg2rad);
                    hRecoPhiLambdaZCut->Fill(lambdaCand.getPhi() * deg2rad);

                    primVertexBeforeDecayVertex++;
                }
                else
                {
                    decayVertexBeforePrimVertex++;
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
    h021->Write();
    h023->Write();
    h03->Write();
    h031->Write();
    h033->Write();
    h04->Write();
    h040->Write();
    h041->Write();
    h044->Write();
    h05->Write();
    h050->Write();
    h051->Write();
    h054->Write();
    h08->Write();
    h081->Write();
    h08->Write();

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
