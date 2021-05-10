#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuildeer(const std::vector<HParticleCand> &particleCands) : fParticleCands(particleCands),
                                                                                 fVerbose(0)
{
    // Loop over all input tracks
    for (int numInputCands = 0, numInputCands < fParticleCands.size(), numInputCands++)
    {

        HParticleCand *inputParticleCand = particleCands[numInputCands];
        estimateCovarianceMatrix(inputParticleCand);
    }

    // Perform the analysis
}

void HDecayBuilder::estimateCovarianceMatrix(HParticleCand *cand)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HDecayBuildeer::estimateCovarianceMatrix() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    if (fMomDepErrors == true)
    {
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
        if (cand->getGeantPID() == 14) //Proton found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                               RErrP->Eval(mom), ZErrP->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, candidate, errors, 938.272);
            fProtons.push_back(candidate);
        }
        else if (cand->getGeantPID() == 9) // Pion found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                               RErrPi->Eval(mom), ZErrPi->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, candidate, errors, 139.570);
            fPions.push_back(candidate);
        }
        else if (cand->getGeantPID() == 11) // Kaon found
        {
            Double_t mom = cand->getGeantTotalMom();
            double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                               RErrK->Eval(mom), ZErrK->Eval(mom)}; //momentum dependent error estimates
            FillData(cand, candidate, errors, 493.7);
            fKaons.push_back(candidate);
        }
        else
            continue;
    }

    if (fFixedErrors == true)
    {
        if (cand->getGeantPID() == 14) //Proton found
        {
            //std::cout << "Parent ID: " << cand->getGeantParentPID() << std::endl;
            if (selectHadrons(cand) == true)
            {
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                fProtons.push_back(candidate);
            }
        }
        else if (cand->getGeantPID() == 9) // Pion found
        {
            if (selectHadrons(cand) == true)
            {
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                fPions.push_back(candidate);
            }
        }
        else if (cand->getGeantPID() == 11) // Kaon found
        {
            if (selectHadrons(cand) == true)
            {
                double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                                   1.404, 2.723};
                FillData(cand, candidate, errors, 493.7);
                fKaons.push_back(candidate);
            }
        }
        else
            continue;
    }
}

void HDecayBuildeer::FillData(particleCand cand, HRefitCand &outcand, double arr[], double mass)
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