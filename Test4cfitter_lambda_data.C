//errors from Waleed's presentation from Aug 28
//no hvirtualcand.h in hydra
/*
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
*/

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "h4cfitter.h"

#include "/lustre/hades/user/jrieger/pp35_data_4charged/forJana/KParticleCand.h"

using namespace std;
//using namespace Particle;

void FillData(KParticleCand* cand, HRefitCand& outcand, double arr[],
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
/*
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
*/

//Int_t Test4cfitter_lambda_data(TString infileList = "/lustre/hades/user/jrieger/pp35_data_4charged/forJana/pp35_skimmed_68*.root", Int_t nEvents = -1)
Int_t Test4cfitter_lambda_data(TString infileList = "/lustre/hades/user/jrieger/pp35_data_4charged/forJana/pp35_skimmed_*.root", Int_t nEvents = -1)
{

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("test4cFit_lambda_data_converged.root", "recreate");

    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 100);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" counts ");
    TH1F *h022 = (TH1F*)h02->Clone("hChi2_probCut");
    h022 -> SetLineColor(kGreen);
    TH1F *h023 = (TH1F*)h02->Clone("hChi2_converged");
    h022 -> SetLineColor(kBlue);

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" counts ");

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" events ");
    h04 -> SetLineColor(kRed);
    TH1F *h042 = (TH1F*)h04->Clone("hLambdaMassPostFit_probCut");
    h042 -> SetLineColor(kGreen);
    TH1F *h043 = (TH1F*)h04->Clone("hLambdaMassPostFit_converged");
    h043 -> SetLineColor(kBlue);

    TH1F* h05 = new TH1F("hPull", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" counts ");
    TH1F *h052 = (TH1F*)h05->Clone("hPull_probCut");
    h052 -> SetLineColor(kGreen);
    TH1F *h053 = (TH1F*)h05->Clone("hPull_converged");
    h053 -> SetLineColor(kBlue);
    
    TH1F* h06 = new TH1F("hTotMomPreFit", "", 100, 3800, 4800);
    h06->SetXTitle(" p [MeV/c]");
    h06->SetYTitle(" events ");
    
    TH1F* h07 = new TH1F("hTotMomPostFit", "", 100, 3800, 4800);
    h07->SetXTitle(" p [MeV/c]");
    h07->SetYTitle(" events ");
    h07 -> SetLineColor(kRed);
    TH1F *h072 = (TH1F*)h07->Clone("hTotMomPostFit_probCut");
    h072 -> SetLineColor(kGreen);
    TH1F *h073 = (TH1F*)h07->Clone("hTotMomPostFit_converged");
    h072 -> SetLineColor(kBlue);

    TH1F* h08 = new TH1F("hNIterations", "", 4, 0, 5);
    h07->SetXTitle(" Iteration");
    h07->SetYTitle(" events ");
    
    
    // -----------------------------------------------------------------------
    
    TChain *t = new TChain("analysis");
    // Add all files here
    t->Add(infileList);
    
    TFile *pid_cuts = new TFile("~/Hades/hydra/pp_pKLambda/pid_cut.root");
    TCutG *cut_proton = (TCutG*)pid_cuts->Get("protons");
    TCutG *cut_pion = (TCutG*)pid_cuts->Get("pions");
    TCutG *cut_kaon = (TCutG*)pid_cuts->Get("kaons");
    
    TClonesArray *cands = new TClonesArray("KParticleCand");
    t->GetBranch("KParticleCand")->SetAutoDelete(kFALSE);
    t->SetBranchAddress("KParticleCand", &cands);

    Long64_t nevts = t->GetEntries();
    if (nEvents <= 0 || nEvents > nevts)
        nEvents = nevts;
    // event loop
    for (Long64_t ev = 0; ev < nEvents; ev++)
    {
        //HTool::printProgress(i, nEvents, 1, "Analysing evt# :");
	
        t->GetEntry(ev);
	
	std::vector<KParticleCand*> protons, kaons, pions;
	std::vector<HRefitCand> protons_fit, pions_fit, kaons_fit;
	
        // track loop
        for (int i = 0; i < cands->GetEntriesFast(); i++)
        {
            KParticleCand *cand = (KParticleCand *)cands->At(i);
	    
	    //do pid
	    Int_t pid = 0;
	    
	    if (cand->getCharge() == 0)
                continue;
		
	    float pq = cand->getMomentum() * cand->getCharge();
	    float beta = cand->getBeta();
	    
	    if( cut_proton->IsInside(pq,beta) ) pid = 14;
	    else if ( cut_pion->IsInside(pq,beta) ) pid = 9;
	    else if ( cut_kaon->IsInside(pq,beta) ) pid = 11;
	    else continue;
	    
	    HRefitCand candidate(cand);
            // proton pdg==14, pion pdg==9
            // error values obtained from resoultion plots
            if (pid == 14)
            {    
	        protons.push_back(cand);
                double errors[] = {2.7 * 1e-5, 4.7 * 1e-3, 1.1 * 1e-2,
                                   2.4, 4.7};
                FillData(cand, candidate, errors, 938.272);
                protons_fit.push_back(candidate);
            }
            else if (pid == 9)
            {
                pions.push_back(cand);
		double errors[] = {9.6 * 1e-5, 1.1 * 1e-2, 2.7 * 1e-2,
                                   6., 11.};
                FillData(cand, candidate, errors, 139.570);
                pions_fit.push_back(candidate);
            }
            else if (pid == 11)
            {
                kaons.push_back(cand);
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
	
	for (size_t k = 0; k < protons_fit.size(); k++)
        {
	    HRefitCand cand1 = protons_fit[k];
	   // cout << " proton1 " << endl;
	for (size_t l = 0; l < kaons_fit.size(); l++)
        {
	    HRefitCand cand2 = kaons_fit[l];
	   // cout << " kaon " << endl;
        for (size_t n = 0; n < protons.size(); n++)
        {
            HRefitCand cand3 = protons_fit[n];
	    //cout << " proton2 " << endl;
            for (size_t m = 0; m < pions_fit.size(); m++)
            {
                HRefitCand cand4 = pions_fit[m];
		//cout << " pion " << endl;
                // mass prefit
		TLorentzVector lambda = cand3 + cand4;
                h01->Fill(lambda.M());
		//cout << " lambda mass "<< lambda.M() << endl;
		
		TLorentzVector all = lambda + cand1 + cand2;
		h06 -> Fill(all.P());
		//cout << " all mom "<< all.P() << endl;
		
		
                // ---------------------------------------------------------------------------------
                // begin kinfit here
                // ---------------------------------------------------------------------------------
                std::vector<HRefitCand> cands;
                cands.push_back(cand1);
                cands.push_back(cand2);
                cands.push_back(cand3);
                cands.push_back(cand4);
		//cout << " cands pushed back " << endl;
		
		TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);
		
		//cout << " ini fitter: " << endl;
                H4cFitter fitter(cands, ppSystem);
		//cout << " ini done " << endl;
                fitter.add4Constraint();
		//cout << " constraint added " << endl;
                if(fitter.fit()){

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
		
		if(fitter.getProb()>0.01){
			h022->Fill(fitter.getChi2());
            h042->Fill(lambda_fit.M());
			h072->Fill(all_fit.P());

            // get Pull example (1/P for the fitted proton)
            h052->Fill(fitter.getPull(0));
		}

        if(fitter.isConverged()){
			h023->Fill(fitter.getChi2());
            h043->Fill(lambda_fit.M());
			h073->Fill(all_fit.P());

            // get Pull example (1/P for the fitted proton)
            h053->Fill(fitter.getPull(0));

            h08->Fill(fitter.getIteration());
		}

		}
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
    h022->Write();
    h023->Write();
    h03->Write();
    h04->Write();
    h042->Write();
    h043->Write();
    h05->Write();
    h052->Write();
    h053->Write();
    h06->Write();
    h07->Write();
    h072->Write();
    h073->Write();
    h08->Write();
    outfile->Close();

    return 0;
} // end of the macro
