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

Int_t analysis(TString infileList = "pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents = 100000)
{

    TStopwatch timer;
    timer.Start();

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("vertexfit_test.root", "recreate");

    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" Counts ");

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" Counts ");

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" Events ");
    
    TH1F* h17 = new TH1F("hLambdaMassPostFitCut", "", 100, 1070, 1170);
    h17->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h17->SetYTitle(" Events ");
    
    TH1F* h05 = new TH1F("hPullPInv", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" Counts ");    
    
    TH1F* h06 = new TH1F("hPullTheta", "", 100, -5, 5);
    h06->SetXTitle("Pull(#theta)");
    h06->SetYTitle(" Counts ");  
    
    TH1F* h07 = new TH1F("hPullPhi", "", 100, -5, 5);
    h07->SetXTitle("Pull(#phi)");
    h07->SetYTitle(" Counts ");  
       
    TH1F* h08 = new TH1F("hPullR", "", 100, -5, 5);
    h08->SetXTitle("Pull(R)");
    h08->SetYTitle(" Counts ");      
    
    TH1F* h09 = new TH1F("hPullZ", "", 100, -5, 5);
    h09->SetXTitle("Pull(Z)");
    h09->SetYTitle(" Counts ");  
    
    TH1F* h10 = new TH1F("hSuccessfulConvergence", "", 100, 1070,1170);
    h10->SetXTitle("M_{p#pi^{-}} [MeV/c^{2}]");
    h10->SetYTitle(" Counts ");    
    
    TH1F* h11 = new TH1F("hIterations", "", 100, 0, 20);
    h11->SetXTitle("Number of Iterations");
    h11->SetYTitle(" Counts ");
    
    // ------------------------- Histos after probability cut    --------------------
    
    TH1F* h12 = new TH1F("hPullPInvCut", "", 100, -5, 5);
    h12->SetXTitle("Pull(1/P_{p})");
    h12->SetYTitle(" Counts ");    
    
    TH1F* h13 = new TH1F("hPullThetaCut", "", 100, -5, 5);
    h13->SetXTitle("Pull(#theta)");
    h13->SetYTitle(" Counts ");  
    
    TH1F* h14 = new TH1F("hPullPhiCut", "", 100, -5, 5);
    h14->SetXTitle("Pull(#phi)");
    h14->SetYTitle(" Counts ");  
       
    TH1F* h15 = new TH1F("hPullRCut", "", 100, -5, 5);
    h15->SetXTitle("Pull(R)");
    h15->SetYTitle(" Counts ");      
    
    TH1F* h16 = new TH1F("hPullZCut", "", 100, -5, 5);
    h16->SetXTitle("Pull(Z)");
    h16->SetYTitle(" Counts ");  

    TH1F* hMinDistanceBetweenCands = new TH1F("hMinDistanceBetweenCands", "", 100,0,100);
    hMinDistanceBetweenCands->SetXTitle("Minimum Distance Between Tracks / cm");
    hMinDistanceBetweenCands->SetYTitle("Occurrence");
    
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
    if (!loop.setInput("-*,+HParticleCandSim,+HGeantKine"))
    {
        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory* catParticle = loop.getCategory("HParticleCandSim");
    if (!catParticle)
    {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }
    HCategory* catGeant = loop.getCategory("HGeantKine");
    if (!catGeant)
    {
        std::cout << "No kineCat in input!" << std::endl;
        exit(1);
    }

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
	

        std::vector<HRefitCand> protons, pions, kaons;
        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim* cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack()) continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed)) continue;

            HRefitCand candidate(cand);
            
	    // select particles based on MC info
            // proton pdg==14, pion pdg==9, k+ pdg==11, lambda pdg==18
            // error values obtained from resoultion plots
	    
            if (cand->getGeantPID() == 14 && cand->getGeantParentPID() == 18) //Proton found
            {
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9 && cand->getGeantParentPID() == 18) // Pion found
            {
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }          
	    else if (cand->getGeantPID() == 11) // Kaon found
            {
                double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                                   1.404, 2.723};
                FillData(cand, candidate, errors, 493.7);
                kaons.push_back(candidate);
            }
            else
                continue;
        } // end track loop


        // Testing the minimum distance between the two tracks
        for (size_t n = 0; n < protons.size(); n++)
        {
            TVector3 base1, base2, dir1, dir2;

            HRefitCand cand1 = protons[n];
            
            base1.SetXYZ(0,0,0); 
            dir1.SetXYZ(0,0,0);
            
            base1.SetXYZ(cand1.getR() * std::cos(cand1.Phi() + TMath::PiOver2()),
                          cand1.getR() * std::sin(cand1.Phi() + TMath::PiOver2()),
                          cand1.getZ());

            dir1.SetXYZ(std::sin(cand1.Theta()) * std::cos(cand1.Phi()),
                        std::sin(cand1.Theta()) * std::sin(cand1.Phi()),
                        std::sin(cand1.Theta()));

            for (size_t m = 0; m < pions.size(); m++)
            {

                HRefitCand cand2 = pions[m];

                base2.SetXYZ(0,0,0);
                dir2.SetXYZ(0,0,0);

                base2.SetXYZ(cand2.getR() * std::cos(cand2.Phi() + TMath::PiOver2()),
                              cand2.getR() * std::sin(cand2.Phi() + TMath::PiOver2()),
                              cand2.getZ());

                dir2.SetXYZ(std::sin(cand2.Theta()) * std::cos(cand2.Phi()),
                            std::sin(cand2.Theta()) * std::sin(cand2.Phi()),
                            std::sin(cand2.Theta()));

                            
                            // Calculate minimum distance between the tracks

                            double d = std::fabs((dir1.Cross(dir2)).Dot((base1 - base2)));

                            hMinDistanceBetweenCands->Fill(d);


			     for (size_t n = 0; n < protons.size(); n++)
        			{
            				HRefitCand cand1 = protons[n];

            }
        }

 /*       // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------
        for (size_t n = 0; n < protons.size(); n++)
        {for (size_t m = 0; m < pions.size(); m++)
            {
                HRefitCand cand2 = pions[m];

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
        for (size_t n = 0; n < protons.size(); n++)
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
        } 
	
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
    hMinDistanceBetweenCands->Write();

    outfile->Close();

    return 0;
} // end of the macro
