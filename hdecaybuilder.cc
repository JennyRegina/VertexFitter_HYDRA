#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuildeer(const std::vector<HParticleCand> &particleCands) : fParticleCands(particleCands),
                                                                              fVerbose(0)
{
    
    for(int numInputCands=0, numInputCands<fParticleCands.size(), numInputCands++){


        HParticleCand inputParticleCand = particleCands[numInputCands];
        FillData(inputParticleCand, HRefitCand & outcand, double arr[], double mass);

    }

}

void HDecayBuildeer::FillData(particleCand, HRefitCand & outcand, double arr[], double mass)
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