//errors from Waleed's presentation from Aug 28

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

#include "h4cfitter.h"

using namespace std;
using namespace Particle;

void FillData(HParticleCand* cand, HRefitCand& outcand, double arr[],
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

Bool_t selectHadrons(HParticleCand* pcand)
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

    if (!test) return kFALSE;

    if (test) test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE;

    return test;
}


Int_t Test4cfitter_lambda(TString infileList = "/lustre/hades/user/jrieger/pp_pKLambda/sim/pp_pKLambda_3500_hgeant1_dst_apr12.root", Int_t nEvents = 50000)
{

    TStopwatch timer;
    timer.Start();

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("test4cFit_lambda.root", "recreate");

    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" counts ");

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" counts ");

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" events ");

    TH1F* h05 = new TH1F("hPull", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" counts ");
    
    TH1F* h06 = new TH1F("hTotMomPreFit", "", 100, 3800, 4800);
    h06->SetXTitle(" p [MeV/c]");
    h06->SetYTitle(" events ");
    
    TH1F* h07 = new TH1F("hTotMomPostFit", "", 100, 3800, 4800);
    h07->SetXTitle(" p [MeV/c]");
    h07->SetYTitle(" events ");
    
    
    // -----------------------------------------------------------------------

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << infileList.Data() << endl;
        return 1;
    }

    // select categories here
    if(!loop.setInput("-*,+HParticleCandSim,+HGeantKine")) {
        cout<<"READBACK: ERROR : cannot read input !"<<endl;
	    exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory *catParticle = loop.getCategory("HParticleCandSim"); 
    if (!catParticle) { std::cout<<"No particleCat in input!"<<std::endl; exit(1);}
    HCategory *catGeant    = loop.getCategory("HGeantKine");
    if (!catGeant) { std::cout<<"No kineCat in input!"<<std::endl; exit(1);}

    Int_t entries = loop.getEntries();
    if (nEvents < entries && nEvents >= 0) entries = nEvents;

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


        std::vector<HParticleCandSim*> protons1, kaons, protons2, pions;
	std::vector<HRefitCand> protons1_fit, protons2_fit, pions_fit, kaons_fit;
        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim* cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack()) continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed)) continue;
	    
	    if (cand->getGeantPID()==14&&cand->getGeantTrack()==1) protons1.push_back(cand);
            else if (cand->getGeantPID()==9&&cand->getGeantParentTrackNum()==3) pions.push_back(cand);
            else if (cand->getGeantPID()==11&&cand->getGeantTrack()==2) kaons.push_back(cand);
            else if (cand->getGeantPID()==14&&cand->getGeantParentTrackNum()==3) protons2.push_back(cand);

            HRefitCand candidate(cand);
            // select particles based on MC info
            // proton pdg==14, pion pdg==9
            // error values obtained from resoultion plots
            if (cand->getGeantPID() == 14&&cand->getGeantParentTrackNum()==3)
            {
                double errors[] = {2.7 * 1e-5, 4.7 * 1e-3, 1.1 * 1e-2,
                                   2.4, 4.7};
                FillData(cand, candidate, errors, 938.272);
                protons2_fit.push_back(candidate);
            }
	    else if (cand->getGeantPID() == 14&&cand->getGeantTrack()==1)
            {
                double errors[] = {2.7 * 1e-5, 4.7 * 1e-3, 1.1 * 1e-2,
                                   2.4, 4.7};
                FillData(cand, candidate, errors, 938.272);
                protons1_fit.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9)
            {
                double errors[] = {9.6 * 1e-5, 1.1 * 1e-2, 2.7 * 1e-2,
                                   6., 11.};
                FillData(cand, candidate, errors, 139.570);
                pions_fit.push_back(candidate);
            }
            else if (cand->getGeantPID() == 11)
            {
                double errors[] = {1.9 * 1e-5, 5.0 * 1e-3, 6.3 * 1e-3,
                                   2.3, 4.7};
                FillData(cand, candidate, errors, 493.577);
                kaons_fit.push_back(candidate);
            }
            else
                continue;
        } // end track loop
	

        // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------
	
	for (size_t k = 0; k < protons1_fit.size(); k++)
        {
	    HRefitCand cand1 = protons1_fit[k];
	   // cout << " proton1 " << endl;
	for (size_t l = 0; l < kaons_fit.size(); l++)
        {
	    HRefitCand cand2 = kaons_fit[l];
	   // cout << " kaon " << endl;
        for (size_t n = 0; n < protons2_fit.size(); n++)
        {
            HRefitCand cand3 = protons2_fit[n];
	    //cout << " proton2 " << endl;
            for (size_t m = 0; m < pions_fit.size(); m++)
            {
                HRefitCand cand4 = pions_fit[m];
		cout << " pion " << endl;
                // mass prefit
		TLorentzVector lambda = cand3 + cand4;
                h01->Fill(lambda.M());
		cout << " lambda mass "<< lambda.M() << endl;
		
		TLorentzVector all = lambda + cand1 + cand2;
		h06 -> Fill(all.P());
		cout << " all mom "<< all.P() << endl;
		
		
                // ---------------------------------------------------------------------------------
                // begin kinfit here
                // ---------------------------------------------------------------------------------
                std::vector<HRefitCand> cands;
                cands.push_back(cand1);
                cands.push_back(cand2);
                cands.push_back(cand3);
                cands.push_back(cand4);
		cout << " cands pushed back " << endl;
		
		TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);
		
		cout << " ini fitter: " << endl;
                H4cFitter fitter(cands, ppSystem);
		cout << " ini done " << endl;
                fitter.add4Constraint();
		cout << " constraint added " << endl;
                fitter.fit();

                // get fitted objects fittedcand1 and fittedcand2
                HRefitCand fcand1 = fitter.getDaughter(0); // proton
                HRefitCand fcand2 = fitter.getDaughter(1); // kaon
                HRefitCand fcand3 = fitter.getDaughter(2); // proton
                HRefitCand fcand4 = fitter.getDaughter(3); // pion

                h02->Fill(fitter.getChi2());
                h03->Fill(fitter.getProb());
		TLorentzVector lambda_fit = fcand3 + fcand4;
                h04->Fill(lambda_fit.M());
		TLorentzVector all_fit = lambda_fit + fcand1 + fcand2;
                h07->Fill(all_fit.P());

                // get Pull example (1/P for the fitted proton)
                h05->Fill(fitter.getPull(0));
		
                // ---------------------------------------------------------------------------------
            }
        }
	}
	}
        // -----------------------------------------------------------------------

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
    outfile->Close();

    return 0;
} // end of the macro
