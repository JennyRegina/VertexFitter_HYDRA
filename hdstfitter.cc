#include "hdstfitter.h"

HDSTFitter::HDSTFitter(TString infileList, bool includeFw = false, bool momDepErrors=false, Int_t nEvents=-1) : fInfilelist(infileList),
                                                                            fIncludeFw(includeFw),
                                                                            fMomDepErrors(momDepErrors)
{
}

void HDSTFitter::estimateCov(Int_t pid, Double_t mom=0, double& cov[])
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HDecayBuildeer::estimateCovarianceMatrix() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    if (fMomDepErrors == true)
    {
        //Momentum dependent uncertainty estimation input for HADES particle candidates
        TFile *momErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepMomErrors_allParticles_recoMom_errFunc.root", "read");
        TF1 *momErrP = (TF1 *)momErr->Get("f_pP");
        TF1 *momErrPi = (TF1 *)momErr->Get("f_pPi");
        TF1 *momErrK = (TF1 *)momErr->Get("f_pK");
        TF1 *momErrEm = (TF1 *)momErr->Get("f_pEm");

        TFile *thtErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepThtErrors_recoMom_errFunc.root", "read");
        TF1 *thtErrP = (TF1 *)thtErr->Get("f_thtP");
        TF1 *thtErrPi = (TF1 *)thtErr->Get("f_thtPi");
        TF1 *thtErrK = (TF1 *)thtErr->Get("f_thtK");
        TF1 *thtErrEm = (TF1 *)thtErr->Get("f_thtEm");

        TFile *phiErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepPhiErrors_recoMom_errFunc.root", "read");
        TF1 *phiErrP = (TF1 *)phiErr->Get("f_phiP");
        TF1 *phiErrPi = (TF1 *)phiErr->Get("f_phiPi");
        TF1 *phiErrK = (TF1 *)phiErr->Get("f_phiK");
        TF1 *phiErrEm = (TF1 *)phiErr->Get("f_phiEm");

        TFile *RZErr_PPi = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_PPi_errFunc.root", "read");
        TFile *RZErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_PK_errFunc.root", "read");
        TFile *RZErr_PEm = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_recoMom_PEm_errFunc.root", "read");
        TF1 *RErrP = (TF1 *)RZErr_PEm->Get("f_RP");
        TF1 *RErrPi = (TF1 *)RZErr_PPi->Get("f_RPi");
        TF1 *RErrK = (TF1 *)RZErr_PK->Get("f_RK");
        TF1 *RErrEm = (TF1 *)RZErr_PEm->Get("f_REm");
        TF1 *ZErrP = (TF1 *)RZErr_PEm->Get("f_ZP");
        TF1 *ZErrPi = (TF1 *)RZErr_PPi->Get("f_ZPi");
        TF1 *ZErrK = (TF1 *)RZErr_PK->Get("f_ZK");
        TF1 *ZErrEm = (TF1 *)RZErr_PEm->Get("f_ZEm");

        // FwDet protons
        TFile *momErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepMomErrors_fw_recoMom_errFunc.root", "read");
        TF1 *momErrP_fw = (TF1 *)momErr_fw->Get("f_pP");
        TFile *thtErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepThtErrors_fw_recoMom_errFunc.root", "read");
        TF1 *thtErrP_fw = (TF1 *)thtErr_fw->Get("f_thtP");
        TFile *phiErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepPhiErrors_fw_recoMom_errFunc.root", "read");
        TF1 *phiErrP_fw = (TF1 *)phiErr_fw->Get("f_phiP");
        TFile *RZErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_fw_recoMom_errFunc.root", "read");
        TF1 *RErrP_fw = (TF1 *)RZErr_fw->Get("f_RP");
        TF1 *ZErrP_fw = (TF1 *)RZErr_fw->Get("f_ZP");

        if (pid == 14) cov[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                               RErrP->Eval(mom), ZErrP->Eval(mom)};
        else if (pid == 9) cov[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                               RErrPi->Eval(mom), ZErrPi->Eval(mom)};
        else if (pid == 11) cov[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                               RErrK->Eval(mom), ZErrK->Eval(mom)};
        else if (pid == 3) cov[] = {momErrEm->Eval(mom), thtErrEm->Eval(mom), phiErrEm->Eval(mom),
                               RErrEm->Eval(mom), ZErrEm->Eval(mom)};
        else cout<<"No momentum dependent error estimate available for this pid"<<endl;
    }

    if (fFixedErrors == true) //which momentum/setup is this for?
    {
        if (pid == 14) cov[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                               1.188, 2.652};
        else if (pid == 9) cov[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                               4.006, 7.629};
        else if (pid == 11) cov[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                               1.404, 2.723};
        else cout<<"No error estimate available for this pid"<<endl;
    }
}

void HDSTFitter::selectCandidates()
{
    // for each event there are a number of tracks
    Int_t ntracks = catParticle->getEntries();
    if(fIncludeFw) Int_t nFwTracks = catFwParticle->getEntries();

    std::vector<HRefitCand> cands_fit[];
    
    for (Int_t j = 0; j < ntracks; j++)
    {
        HParticleCand* cand = HCategoryManager::getObject(cand, catParticle, j);
        // skip ghost tracks (only avalible for MC events)
        if (cand->isGhostTrack()) continue;
        // select "good" tracks
        if (!cand->isFlagBit(Particle::kIsUsed)) continue;

        HRefitCand candidate(cand);

        for(auto it = std::begin(fPids); it != std::end(fPids); ++it) {
            if (cand->getGeantPID()==fPids[it]){
                Double_t mom = cand->P();
                //vector<double> errors;
                //getErrors(fPids[it], mom, errors);
                Double_t errors[];
                etstimateCov(fPid[it], mom, errors);
                FillData(cand, candidate, errors, cand->getMass());
                fCandsFit[it].push_back(candidate);
            } 
        }
    }   // end of HADES track loop

    if(fIncludeFw){ //find correct index of proton HRefitCand vector
        for(Int_t j=0; j<nFwTracks; j++){
            HFwDetCandSim *cand = HCategoryManager::getObject(cand,catFwParticle,j);
            cand->calc4vectorProperties(938.272);
                
            HRefitCand candidate(cand);
            
            // select particles based on MC info
            if (cand->getGeantPID()==14){
                Double_t mom = cand->P();
                double errors[] = {momErrP_fw->Eval(mom), thtErrP_fw->Eval(mom), phiErrP_fw->Eval(mom),
                                RErrP_fw->Eval(mom), ZErrP_fw->Eval(mom)};
                FillDataFw(cand, candidate, errors, 938.272);
                fCandsFit[it that is 14].push_back(candidate);
            }
            else continue;
        }
    } // end fwTrack loop
}

void HDSTFitter::addBuilderTask(TString val, std::vector<Int_t> pids, TLorentzVector lv = (0,0,0,0)){

    fCandsFit.Clear();
    selectCandidates();

    //initialize DecayBuilder

    //Write output category

    } //end of the event loop
}

void HDSTFitter::addFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv = (0,0,0,0), Double_t mm=0){
    
    setPids(pids);
    
    TStopwatch timer;
    timer.Start();

    HLoop loop(kTrue);
    Bool:t ret = loop.addFiles(fInfileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << infileList.Data() << endl;
        return 1;
    }

    // select categories here
    if(fIncludeFw) {
        if(!loop.setInput("-*,+HParticleCandSim,+HFwDetCandSim")) {
            cout<<"READBACK: ERROR : cannot read input !"<<endl;
            exit(1);
        } // read all categories

        HCategory *catFwParticle = loop.getCategory("HFwDetCandSim"); 
        if (!catFwParticle) { std::cout<<"No FWparticleCat in input!"<<std::endl; exit(1);}
    } else {
        if(!loop.setInput("-*,+HParticleCandSim")) {
            cout<<"READBACK: ERROR : cannot read input !"<<endl;
            exit(1);
        } // read all categories
    }

    loop.printCategories();
    loop.printChain();

    HCategory *catParticle = loop.getCategory("HParticleCandSim"); 
    if (!catParticle) { std::cout<<"No particleCat in input!"<<std::endl; exit(1);}

    Int_t entries = loop.getEntries();
    if(nEvents > entries || nEvents <= 0 ) nEvents = entries;

    cout<<"events: "<<nEvents<<endl;

    // start of the event loop
    for(Int_t i=1; i<nEvents; i++){
        //----------break if last event is reached-------------
        if(loop.nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
        HTool::printProgress(i,nEvents,1,"Analysing evt# :");

        fCandsFit.Clear();
        selectCandidates();

        //initialize DecayBuilder
        HDecayBuilder builder(fCandsFit, task, fPids, lv, mm);
        builder.buildDecay();

        //Get output particles
    }// end of event loop

    //Write output category
}



void HDSTFitter::FillData(HParticleCand* cand, HRefitCand *outcand, double arr[], double mass)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HDSTFitter::FillData() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand->SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand->setR(cand->getR());
    outcand->setZ(cand->getZ());
    outcand->setCovariance(cov);   
}

void HDSTFitter::FillDataFw(HFwDetCandSim* cand, HRefitCand& outcand, double arr[], double mass)
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);
    
    cand -> calc4vectorProperties(mass);
    
    outcand.SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand.setR(cand->getR());
    outcand.setZ(cand->getZ());
    outcand.setIsForward(true);
    outcand.setCovariance( cov );
}
