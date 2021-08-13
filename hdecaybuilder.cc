#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuilder(std::vector< std::vector<HRefitCand> > &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv, HRefitCand mother, Double_t mass) : fCands(cands), 
                                                                                                                                                            fTask(task),
                                                                                                                                                            fPids(pids),
                                                                                                                                                            fCombiCounter(0),
                                                                                                                                                            fProb(0),
                                                                                                                                                            fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder() -----------------" << std::endl;
    }
    
    setIniSys(lv);
    setMother(mother);
    setMass(mass);
    
    fTotalCombos = 1;
    particleCounter.clear();
    for (size_t i=0; i<fPids.size(); i++){
		fTotalCombos *= fCands[i].size();
		particleCounter.push_back(0);
	}
}

void HDecayBuilder::buildDecay(){
    
    while(fCombiCounter<fTotalCombos){
        cout<<"fill fit cands"<<endl;
        fillFitCands();
        if(fTask = "createNeutral"){
            createNeutralCandidate();
        }else if(fTask = "3C"){
            do3cFit();
        }else if(fTask = "4C"){
			cout<<"4C task received"<<endl;
            do4cFit();
        }else if(fTask = "missMom"){
            doMissMomFit();
        }else{
            cout<<"Task not available"<<endl;
        }
    }
    

}

void HDecayBuilder::createNeutralCandidate()
{  }
/*  
    
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createNeutralCandidate() -----------------" << std::endl;
    }

    bool bestPrimVertexFound=false;
    bool bestDecayVertexFound=false;

    TVector3 decayVertex;
    TVector3 primVertex;
    
    int indexPrimaryProton=-1, indexSecondBestPrimaryProton=-1;
    int indexDecayProton=-1, indexSecondBestDecayProton=-1;

    std::vector<HRefitCand> cands3c;
    cands3c.clear();
    fOutputCands.clear();

    double probPrim = -99999, probSec = -99999;
    double probSecondBestPrim = -99999, probSecondBestDecay = -99999;
    double probDecayVertex_Temp = -1;
    double probPrimVertex_Temp = -1;
    
    // Perform the analysis
    for (size_t n = 0; n < fProtons.size(); n++)
    {
        HRefitCand* cand1Ptr = fProtons[n];
        HRefitCand cand1= *cand1Ptr;
        
        for (size_t m = 0; m < fPions.size(); m++)
        {
            HRefitCand* cand2Ptr = fPions[m]; 
            HRefitCand cand2= *cand2Ptr;               
            
            std::vector<HRefitCand> candsSec;
            candsSec.clear();
            candsSec.push_back(cand1);
            candsSec.push_back(cand2);
            HVertexFinder *vtxFinderSec = new HVertexFinder();
            decayVertex = vtxFinderSec->findVertex(candsSec);
            
            // Kinematic fitting with vertex constraint
            HKinFitter vtxFitterSecCands(candsSec); 

            vtxFitterSecCands.addVertexConstraint();
            vtxFitterSecCands.fit();
            probSecondBestDecay=probSec;
            probSec = vtxFitterSecCands.getProb();

            if (probSec > probDecayVertex_Temp)
            {
                bestDecayVertexFound = true;
                probDecayVertex_Temp = probSec;
                indexSecondBestDecayProton = indexDecayProton;
                indexDecayProton = n;

                cands3c.clear();
                // Get the proton daughter and pass it to the 3C fit later
                cands3c.push_back(vtxFitterSecCands.getDaughter(0));
                // Get the pion daughter and pass it to the 3C fit later
                cands3c.push_back(vtxFitterSecCands.getDaughter(1));
                //std::cout << "vertex position: " << decayVertex.X() << " " << decayVertex.Y() << " " << decayVertex.Z() << std::endl;
            }
            
        }

        for (size_t p = 0; p < fKaons.size(); p++)
        {
            HRefitCand* cand3Ptr = fKaons[p];
            HRefitCand cand3= *cand3Ptr;

            std::vector<HRefitCand> candsPrim;
            candsPrim.clear();
            candsPrim.push_back(cand1);
            candsPrim.push_back(cand3);
            HVertexFinder *vtxFinderPrim = new HVertexFinder();
            primVertex = vtxFinderPrim->findVertex(candsPrim);            
            
            HKinFitter vtxFitterPrimCands(candsPrim);
            vtxFitterPrimCands.addVertexConstraint();
            vtxFitterPrimCands.fit();
            probSecondBestPrim=probPrim;
            probPrim = vtxFitterPrimCands.getProb();

            if (probPrim > probPrimVertex_Temp)
            {
                bestPrimVertexFound = true;
                probPrimVertex_Temp = probPrim;
                indexSecondBestPrimaryProton = indexPrimaryProton;
                indexPrimaryProton = n;                    
                
                fCandsMissPos.push_back(cand1);
                fCandsMissPos.push_back(cand3);
            }
        }
    }

    if (bestPrimVertexFound == true && bestDecayVertexFound == true)
    {

        if (indexDecayProton != indexPrimaryProton && probPrimVertex_Temp>fPrimVertexProbabilityCut && probDecayVertex_Temp>fDecayVertexProbabilityCut)
        {
            HNeutralCandFinder lambdaCandFinder(cands3c);

            lambdaCandFinder.setUsePrimaryVertexInNeutralMotherCalculation(true);

            lambdaCandFinder.setNeutralMotherCandFromPrimaryVtxInfo(primVertex, decayVertex);

            HVirtualCand lambdaCand = lambdaCandFinder.getNeutralMotherCandidate();

            HRefitCand lambdaCandRefit(&lambdaCand);
            lambdaCandRefit.SetXYZM(lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                        std::cos(lambdaCand.getPhi() * deg2rad),
                                    lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                        std::sin(lambdaCand.getPhi() * deg2rad),
                                    lambdaCand.getMomentum() * std::cos(lambdaCand.getTheta() * deg2rad),
                                    1115.683);

            TMatrixD lambdaCov(5, 5);
            lambdaCov = lambdaCandFinder.getCovarianceMatrixNeutralMother();
            lambdaCandRefit.setCovariance(lambdaCov);
            lambdaCandRefit.setR(lambdaCand.getR());
            lambdaCandRefit.setZ(lambdaCand.getZ());
            lambdaCandRefit.SetTheta(lambdaCand.getTheta() * deg2rad);
            lambdaCandRefit.SetPhi(lambdaCand.getPhi() * deg2rad);

            HKinFitter Fitter3c(cands3c, lambdaCandRefit);
            Fitter3c.add3Constraint();
            Fitter3c.setNumberOfIterations(20);

            Fitter3c.fit();

            HRefitCand cand13C = Fitter3c.getDaughter(0); // proton
            createOutputParticle(cand13C);
            HRefitCand cand23C = Fitter3c.getDaughter(1); // pion
            createOutputParticle(cand23C);
            
            HRefitCand lambdaCand3C = Fitter3c.getMother();
            fCandsMissPos.push_back(lambdaCand3C);

        }
    }
}*/

void HDecayBuilder::do4cFit()
{
    HKinFitter Fitter(fFitCands, fIniSys);
    Fitter.add4Constraint();
    cout<<"constraint added"<<endl;
    if(Fitter.fit() && Fitter.getProb()>fProb){
		cout<<"fit successful"<<endl;
        fOutputCands.clear();
        Fitter.getDaughters(fOutputCands);
        /*
        for(iterator it = fPids.begin(); it != fPids.end(); ++it){
            fOutputCands.push_back(Fitter.getDaughter(it));
        }*/
        //fillHistograms();
    }
}

void HDecayBuilder::do3cFit()
{
    HKinFitter Fitter(fFitCands, fMother);
    Fitter.add3Constraint();
    Fitter.fit();
}

void HDecayBuilder::doMissMomFit()
{
    HKinFitter Fitter(fFitCands, fIniSys, fMass);
    Fitter.addMomConstraint();
    Fitter.fit();
}

void HDecayBuilder::fillFitCands()
{
    fFitCands.clear();
    bool doubleParticle = false;
    for (size_t i=0; i<fPids.size(); i++){
		doubleParticle = checkDoubleParticle(i);
        fFitCands.push_back(fCands[i][particleCounter[i]]);
	}
	Int_t a = fPids.size();
	while( particleCounter[a]==fCands[a].size() ){ 
		particleCounter[a] = 0;
		a--;
		if(a<0){
			if(!(fCombiCounter==fTotalCombos)) cout<<"counted wrong: "<<fCombiCounter<<" != "<<fTotalCombos<<endl;
			break;
		}
	}
	particleCounter[a]++;
	fCombiCounter++;
	
    if(doubleParticle) fillFitCands();	//If some particle has been filled more than once into fFitCands, repeat the procedure with the next combination
}

bool checkDoubleParticle(uint i)
{
	for (uint j=0; j<i; i++){
		if((fPids[j]==fPids[i]) && (particleCounter[j]==particleCounter[i])) return true;
	}
	return false;
}
/*
void HDecayBuilder::fillFitCands()
{
    fFitCands.clear();
    //for (size_t i=0; i<fPids.size(); i++){
    for (Int_t i=0; i<fPids.size(); i++){
        fFitCands.push_back(fCands[i][combicounter[fPids[i]]]);
        combicounter[fPids[i]]++;
        if(combicounter[fPids[i]]+1 == fCands[i].size()){
            bool counted = true;
            while(i>=0 && !counted){
                if(combicounter[i]<fCands[i].size()) combicounter[i]++;
                else combicounter[i]=0; i--;
            }
            if(i<0) combicounter[0]=-1;
        }
    }
}

void HDecayBuilder::fillFitCands()
{
    fFitCands.clear();
    //for (size_t i=0; i<fPids.size(); i++){
    for (Int_t i=0; i<fPids.size(); i++){
        fFitCands.push_back(fCands[i][combicounter[fPids[i]]]);
        combicounter[fPids[i]]++;
        if(combicounter[fPids[i]]+1 == fCands[i].size()){
            bool counted = true;
            while(i>=0 && !counted){
                if(combicounter[i]<fCands[i].size()) combicounter[i]++;
                else combicounter[i]=0; i--;
            }
            if(i<0) combicounter[0]=-1;
        }
    }
}*/

/*
void HDecayBuilder::fillFitCands()
{
    fFitCands.clear();
    for (size_t i=0; i<fPids.size(); i++){
        fFitCands.push_back(fCands[i][combicounter[i]]);
        if(i+1 == fPids.size()){
            bool counted = false;
            while(i>=0 && !counted){
                if(combicounter[i]<fCands[i].size()) combicounter[i]++;
                else combicounter[i]=0; i--;
            }
            if(i<0) combicounter[0]=-1;
        }
    }
}*/

/*
void HDecayBuilder::estimateCovarianceMatrix(HParticleCandSim * cand, HRefitCand *refitCand)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HDecayBuildeer::estimateCovarianceMatrix() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    // If both the momentum dependent and fixed errors have been chosen
    // the momentum dependent errors are used since this is more general
    // and needs to be used for the electron fitting 
    
    if(fMomDepErrors == true && fFixedErrors == true){

        fFixedErrors = false;
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

        if (cand->getGeantPID() == 14) //Proton found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                               RErrP->Eval(mom), ZErrP->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 938.272);
            fProtons.push_back(refitCand);
        }
        else if (cand->getGeantPID() == 9) // Pion found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                               RErrPi->Eval(mom), ZErrPi->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 139.570);
            fPions.push_back(refitCand);
        }
        else if (cand->getGeantPID() == 11) // Kaon found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                               RErrK->Eval(mom), ZErrK->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 493.7);
            fKaons.push_back(refitCand);
        }
        else if (cand->getGeantPID() == 3)
        {
            Double_t mom = cand->P();
            double errors[] = {momErrEm->Eval(mom), thtErrEm->Eval(mom), phiErrEm->Eval(mom),
                               RErrEm->Eval(mom), ZErrEm->Eval(mom)};
            FillData(cand, refitCand, errors, 0.51099892);
            fElectrons.push_back(refitCand);
        }
    }

    if (fFixedErrors == true)
    {
        if (cand->getGeantPID() == 14) //Proton found
        {
            double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                               1.188, 2.652};
            FillData(cand, refitCand, errors, 938.272);
            fProtons.push_back(refitCand);
        }
        else if (cand->getGeantPID() == 9) // Pion found
        {
            double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                               4.006, 7.629};
            FillData(cand, refitCand, errors, 139.570);
            fPions.push_back(refitCand);
        }
        else if (cand->getGeantPID() == 11) // Kaon found
        {
            double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                               1.404, 2.723};
            FillData(cand, refitCand, errors, 493.7);
            fKaons.push_back(refitCand);
        }
    }
}*/

void HDecayBuilder::createOutputParticle(HRefitCand refitCand)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputParticle() -----------------" << std::endl;
    }
    HParticleCand *newParticle;

    newParticle->setPhi(refitCand.Theta());//!!!
    newParticle->setR(refitCand.getR());
    newParticle->setZ(refitCand.getZ());
    newParticle->setMomentum(refitCand.P());

    fOutputCands.push_back(newParticle);
}

void HDecayBuilder::createOutputCategory()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputCategory() -----------------" << std::endl;
    }

    HCategoryManager catManager;
    HCategory *cat = catManager.addCategory(1, "HParticleCandSimAfterFit");
}
/*
void HDecayBuilder::fillHistograms()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::writeHistograms() -----------------" << std::endl;
    }

    hmLam_prefit->Fill((fFitCands[2]+fFitCands[3]).M());
    hmLam_post4C->Fill((fOutputCands[2]+fOutputCands[3]).M());
}*/
