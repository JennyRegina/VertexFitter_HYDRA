#include "hdstfitter.h"

HDSTFitter::HDSTFitter(bool includeFw = false, bool momDepErrors=false) : fPIncludeFw(includeFw),
                                                                            fMomDepErrors(momDepErrors),
                                                                            fVerbose(0)
{
}


void HDSTFitter::selectCandidates()
{
    TStopwatch timer;
    timer.Start();

    HLoop loop(kTrue);
    Bool:t ret = loop.addFiles(infileList);
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

            for(auto it = std::begin(pids); it != std::end(pids); ++it) {
                if (cand->getGeantPID()==pids[it]){
                    Double_t mom = cand->P();
                    vector<double> errors;
                    getErrors(pids[it], mom, errors);
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
                    protons_fw.push_back(cand);
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
    
    fCandsFit.Clear();
    selectCandidates();

    //initialize DecayBuilder
    HDecayBuilder builder(fCandsFit, task, pids, lv, mm);
    builder.buildDecay();

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
