#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuilder(std::vector<KParticleCand *> particleCands, std::vector<int> pids) : fParticleCands(particleCands),
                                                                              fVerbose(0)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder() -----------------" << std::endl;
    }

    fProtons.clear();
    fPions.clear();
    fKaons.clear();
    fPosPions.clear();
    fNegKaons.clear();
    fElectrons.clear();
    fPositrons.clear();

    fCandsMissPos.clear();
    
    // Loop over all input tracks
    for (int numInputCands = 0; numInputCands < (int)fParticleCands.size(); numInputCands++)
    {

        KParticleCand *inputParticleCand = particleCands[numInputCands];
        HRefitCand *refitCand = new HRefitCand(inputParticleCand);
        int myPid = pids[numInputCands];

        estimateCovarianceMatrix(inputParticleCand, refitCand, myPid);
    }
    
    // Performs the loop over the particles in the event and findx vertices and the fit
    createNeutralCandidate();

    // Finds electrons and fits them
    if(fFitElectrons == true){
        
        fitElectrons();
    
    }
}

void HDecayBuilder::createNeutralCandidate(){    
    
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

            KParticleCand lambdaCand = lambdaCandFinder.getNeutralMotherCandidate();

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
}

void HDecayBuilder::fitElectrons()
{

    for (size_t o = 0; o < fElectrons.size(); o++)
    {
        HRefitCand* cand5Ptr = fElectrons[o];
        HRefitCand cand5 = *cand5Ptr;

        std::vector<HRefitCand> cands = fCandsMissPos;
        cands.push_back(cand5);

        // TODO: adjust so the user can set the ppSystem
        TLorentzVector ppSystem(0, 0, 5356.7199, 2 * 938.272 + 4500); //initial system
        HKinFitter fitter(cands, ppSystem, 0.51099892);
        //fitter.setVerbosity(verb);
        fitter.addMomConstraint();

        if (fitter.fit() == true)
        {

            // get fitted objects fittedcand1 and fittedcand2
            HRefitCand fcand1 = fitter.getDaughter(0);           // proton
            HRefitCand fcand2 = fitter.getDaughter(1);           // kaon
            HRefitCand fcand3 = fitter.getDaughter(2);           // lambda
            HRefitCand fcand4 = fitter.getDaughter(3);           // electron
            TLorentzVector fcand5 = fitter.getMissingDaughter(); // positron

            //sigma mass after fits
            TLorentzVector sigma = fcand3 + fcand4 + fcand5;
        }

    } //electron loop
}

void HDecayBuilder::estimateCovarianceMatrix(KParticleCand * cand, HRefitCand *refitCand, int myPid)
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

        if (myPid == 14) //Proton found
        {
            Double_t mom = cand->getMomentum();
            double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                               RErrP->Eval(mom), ZErrP->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 938.272);
            fProtons.push_back(refitCand);
        }
        else if (myPid == 9) // Pion found
        {
            Double_t mom = cand->getMomentum();
            double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                               RErrPi->Eval(mom), ZErrPi->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 139.570);
            fPions.push_back(refitCand);
        }
        else if (myPid == 11) // Kaon found
        {
            Double_t mom = cand->getMomentum();
            double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                               RErrK->Eval(mom), ZErrK->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, refitCand, errors, 493.7);
            fKaons.push_back(refitCand);
        }
        else if (myPid == 3)
        {
            Double_t mom = cand->getMomentum();
            double errors[] = {momErrEm->Eval(mom), thtErrEm->Eval(mom), phiErrEm->Eval(mom),
                               RErrEm->Eval(mom), ZErrEm->Eval(mom)};
            FillData(cand, refitCand, errors, 0.51099892);
            fElectrons.push_back(refitCand);
        }
    }

    if (fFixedErrors == true)
    {
        if (myPid == 14) //Proton found
        {
            double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                               1.188, 2.652};
            FillData(cand, refitCand, errors, 938.272);
            fProtons.push_back(refitCand);
        }
        else if (myPid == 9) // Pion found
        {
            double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                               4.006, 7.629};
            FillData(cand, refitCand, errors, 139.570);
            fPions.push_back(refitCand);
        }
        else if (myPid == 11) // Kaon found
        {
            double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                               1.404, 2.723};
            FillData(cand, refitCand, errors, 493.7);
            fKaons.push_back(refitCand);
        }
    }
}

void HDecayBuilder::FillData(KParticleCand* cand, HRefitCand *outcand, double arr[], double mass)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HDecayBuildeer::FillData() -----------" << std::endl;
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

void HDecayBuilder::createOutputParticle(HRefitCand refitCand)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputParticle() -----------------" << std::endl;
    }
    KParticleCand *newParticle;

    newParticle->setPhi(refitCand.Theta());
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