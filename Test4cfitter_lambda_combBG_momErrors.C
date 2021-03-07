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


Int_t Test4cfitter_lambda_combBG_momErrors(TString infileList = "/lustre/hades/user/jrieger/pp_pKLambda/sim/pp_pKLambda_3500_hgeant1_dst_apr12.root,/lustre/hades/user/jregina/DstTest/pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents = -1)
{

    TStopwatch timer;
    timer.Start();

    //Uncertainty estimation input
    TFile *momErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepMomErrors_allParticles_errFunc.root", "read");
    TF1 *momErrP = (TF1*) momErr->Get("f_pP");
    TF1 *momErrPi = (TF1*) momErr->Get("f_pPi");
    TF1 *momErrK = (TF1*) momErr->Get("f_pK");

    TFile *thtErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_allParticles_errFunc.root", "read");
    TFile *thtErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_PK_errFunc.root", "read");
    TF1 *thtErrP = (TF1*) thtErr->Get("f_thtP");
    TF1 *thtErrPi = (TF1*) thtErr->Get("f_thtPi");
    TF1 *thtErrK = (TF1*) thtErr_PK->Get("f_thtK");

    TFile *phiErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PPi_errFunc.root", "read");    //if 4500 in name this is just a name issue
    TFile *phiErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PK_errFunc.root", "read");
    TF1 *phiErrP = (TF1*) phiErr->Get("f_phiP");
    TF1 *phiErrPi = (TF1*) phiErr->Get("f_phiPi");
    TF1 *phiErrK = (TF1*) phiErr_PK->Get("f_phiK");

    TFile *RZErr_PPi = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PPi_errFunc.root", "read");
    TFile *RZErr_PK= new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PK_errFunc.root", "read");
    TF1 *RErrP = (TF1*) RZErr_PPi->Get("f_RP");
    TF1 *RErrPi = (TF1*) RZErr_PPi->Get("f_RPi");
    TF1 *RErrK = (TF1*) RZErr_PK->Get("f_RK");
    TF1 *ZErrP = (TF1*) RZErr_PPi->Get("f_ZP");
    TF1 *ZErrPi = (TF1*) RZErr_PPi->Get("f_ZPi");
    TF1 *ZErrK = (TF1*) RZErr_PK->Get("f_ZK");

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("test4cFit_lambda_combBG_lr1_it10_momDepErrors_basicChi2Formula_noupdatef.root", "recreate");
    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 100);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" counts ");
    TH1F *h022 = (TH1F*)h02->Clone("hChi2_probCut");
    h022 -> SetLineColor(kGreen);
    TH1F *h023 = (TH1F*)h02->Clone("hChi2_converged");
    h023 -> SetLineColor(kBlue);

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" counts ");
    TH1F *h033 = (TH1F*)h03->Clone("hPChi2_converged");
    h033 -> SetLineColor(kBlue);

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" events ");
    h04 -> SetLineColor(kRed);
    TH1F *h042 = (TH1F*)h04->Clone("hLambdaMassPostFit_probCut");
    h042 -> SetLineColor(kGreen);
    TH1F *h043 = (TH1F*)h04->Clone("hLambdaMassPostFit_converged");
    h043 -> SetLineColor(kBlue);
    TH1F *h044 = (TH1F*)h04->Clone("hLambdaMassPostFit_converged_probCut");
    h043 -> SetLineColor(kCyan);

    TH1F* h05 = new TH1F("hPull", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" counts ");
    TH1F *h052 = (TH1F*)h05->Clone("hPull_probCut");
    h052 -> SetLineColor(kGreen);
    TH1F *h053 = (TH1F*)h05->Clone("hPull_converged");
    h053 -> SetLineColor(kBlue);
    TH1F *h054 = (TH1F*)h05->Clone("hPull_converged_probCut");
    h054 -> SetLineColor(kCyan);
    std::vector<TH1F*> hpulls, hpulls_conv, hpulls_conv_probCut;
    for(int i=0; i<15; i++){
        TH1F *hpull = (TH1F*)h05->Clone("hpull");
        hpulls.push_back(hpull);
        TH1F *hpull_conv = (TH1F*)h053->Clone("hpull_conv");
        hpulls_conv.push_back(hpull_conv);
        TH1F *hpull_conv_probCut = (TH1F*)h05->Clone("hpull_conv_probCut");
        hpulls_conv_probCut.push_back(hpull_conv_probCut);
    }
    
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
    h073 -> SetLineColor(kBlue);
    TH1F *h074 = (TH1F*)h07->Clone("hTotMomPostFit_converged_probCut");
    h074 -> SetLineColor(kCyan);

    TH1F* h08 = new TH1F("hNIterations", "", 10, 0, 10);
    h08->SetXTitle(" Iteration");
    h08->SetYTitle(" events ");
    TH1F *h084 = (TH1F*)h08->Clone("hNIterations_probCut");
    h084 -> SetLineColor(kCyan);
            
    TCanvas *c=new TCanvas("c", "c", 20,20,800,800);
    
    // -----------------------------------------------------------------------

    HLoop loop(kTRUE);
    Bool_t ret = loop.addMultFiles(infileList);
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
    if (nEvents > entries || nEvents < 0) nEvents = entries;
    cout<<"nEvents set to "<<nEvents<<endl;

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


        std::vector<HParticleCandSim*> protons, kaons, pions;
	std::vector<HRefitCand> protons_fit, pions_fit, kaons_fit;
        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim* cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack()) continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed)) continue;
	    
	    if (cand->getGeantPID()==14) protons.push_back(cand);
            else if (cand->getGeantPID()==9) pions.push_back(cand);
            else if (cand->getGeantPID()==11) kaons.push_back(cand);

            HRefitCand candidate(cand);
            // select particles based on MC info
            // proton pdg==14, pion pdg==9
            // error values obtained from resoultion plots
            
            if (cand->getGeantPID() == 14)
            {
                Double_t mom = cand->getGeantTotalMom();
                double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                                   RErrP->Eval(mom), ZErrP->Eval(mom)};
                FillData(cand, candidate, errors, 938.272);
                protons_fit.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9)
            {
                Double_t mom = cand->getGeantTotalMom();
                double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                                   RErrPi->Eval(mom), ZErrPi->Eval(mom)};
                FillData(cand, candidate, errors, 139.570);
                pions_fit.push_back(candidate);
            }
            else if (cand->getGeantPID() == 11)
            {
                Double_t mom = cand->getGeantTotalMom();
                double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                                   RErrK->Eval(mom), ZErrK->Eval(mom)};
                FillData(cand, candidate, errors, 493.577);
                kaons_fit.push_back(candidate);
            }
            else
                continue;
/*
            if (cand->getGeantPID()==14){
	    	//double errors[] = {2.7 * 1e-5, 4.7 * 1e-3, 1.1 * 1e-2, 2.4, 4.7};    //Waleed
		double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3, 1.188, 2.652};    //Rafal
            for(int j = 0; j < 5; j++) errors[j] *= 0.5;
                FillData(cand, candidate, errors, 938.272);
                protons_fit.push_back(candidate);
	    }
            else if (cand->getGeantPID()==9){ 
		//double errors[] = {9.6 * 1e-5, 1.1 * 1e-2, 2.7 * 1e-2, 6., 11.};    //Waleed
		double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2, 4.006, 7.629};    //Rafal
            for(int j = 0; j < 5; j++) errors[j] *= 0.5;
                FillData(cand, candidate, errors, 139.570);
                pions_fit.push_back(candidate);
            }
            else if (cand->getGeantPID()==11){ 
	    	//double errors[] = {3. * 1e-5, 5. * 1e-3, 9. * 1e-3, 2., 5.};    //Guess
            double errors[]  = {1.9 * 1e-5, 5.0 * 1e-3, 6.3 * 1e-3, 2.3, 4.7}; //Waleed
            for(int j = 0; j < 5; j++) errors[j] *= 0.5;
                FillData(cand, candidate, errors, 493.577);
                kaons_fit.push_back(candidate);
	        }
	    else continue;*/
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
                if(fitter.fit(1., 10)){

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

            for(int i=0; i<hpulls.size(); i++){
                TH1F *hpull = hpulls[i];
                hpull->Fill(fitter.getPull(i));
            }

            if(fitter.isConverged()){
                h044->Fill(lambda_fit.M());
			    h074->Fill(all_fit.P());

                // get Pull example (1/P for the fitted proton)
                h054->Fill(fitter.getPull(0));// get Pull
                for(int i=0; i<hpulls_conv_probCut.size(); i++){
                  TH1F *hpull_conv_probCut = hpulls_conv_probCut[i];
                  hpull_conv_probCut->Fill(fitter.getPull(i));
                }

                h084->Fill(fitter.getIteration());
		    }
		}
        if(fitter.isConverged()){
			h023->Fill(fitter.getChi2());
            h033->Fill(fitter.getProb());
            h043->Fill(lambda_fit.M());
			h073->Fill(all_fit.P());

            // get Pull example (1/P for the fitted proton)
            h053->Fill(fitter.getPull(0));

            h08->Fill(fitter.getIteration());

            // get Pull
            for(int i=0; i<hpulls_conv.size(); i++){
                TH1F *hpull_conv = hpulls_conv[i];
                hpull_conv->Fill(fitter.getPull(i));
            }

		}

		}
                // ---------------------------------------------------------------------------------
            }
        }
	}
	}
        // -----------------------------------------------------------------------

    } // end of the events loop

    c->cd();
    for(int i=0; i<hpulls.size(); i++){
        TH1F *hpull = hpulls[i];
        hpull->Draw();
        TString image = Form("pull4cFitter/pull%i.png",i);
        c -> Print(image);
    }

    for(int i=0; i<hpulls_conv.size(); i++){
        TH1F *hpull_conv = hpulls_conv[i];
        hpull_conv->Draw();
        TString image = Form("pull4cFitter_conv/pull%i.png",i);
        c -> Print(image);
    }

    for(int i=0; i<hpulls_conv_probCut.size(); i++){
        TH1F *hpull_conv_probCut = hpulls_conv_probCut[i];
        hpull_conv_probCut->Draw();
        TString image = Form("pull4cFitter_conv_probCut/pull%i.png",i);
        c -> Print(image);
    }

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
    h044->Write();
    h05->Write();
    h052->Write();
    h053->Write();
    h054->Write();
    h06->Write();
    h07->Write();
    h072->Write();
    h073->Write();
    h074->Write();
    h08->Write();
    h084->Write();
    outfile->Close();

    return 0;
} // end of the macro
