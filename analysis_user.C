#include "hades.h"
#include "hloop.h"
#include "htool.h"
#include "hcategorymanager.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hphysicsconstants.h"
#include "hhistmap.h"
#include "hparticletracksorter.h"
#include "henergylosscorrpar.h"


#include "hcategory.h"
#include "hlinearcategory.h"
#include "hrichhit.h"
#include "hrichhitsim.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlepair.h"
#include "hparticlegeantpair.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"

#include "hfwdetcandsim.h"

#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"
#include "hparticlecutrange.h"

#include "TTree.h"

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>
#include <math.h>
/*
#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"
*/
#include "hdstfitter.h"

using namespace std;
using namespace Particle;

Int_t analysis_user(TString infileList="/lustre/hades/user/jrieger/pp_pKLambda/sim/pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents=-1){

    TFile *outfile = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/test_userfit.root","recreate");

    TH1F* hmLam_prefit = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    hmLam_prefit->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    hmLam_prefit->SetYTitle(" events ");
    hmLam_prefit->GetXaxis()->SetTitleSize(0.05);
    hmLam_prefit->GetXaxis()->SetLabelSize(0.05);
    hmLam_prefit->GetYaxis()->SetTitleSize(0.05);
    hmLam_prefit->GetYaxis()->SetLabelSize(0.05);
    hmLam_prefit->SetLineColor(kBlack);
    TH1F *hmLam_post4C = (TH1F*)hmLam_prefit->Clone("hmLam_post4C");
    hmLam_post4C->SetLineColor(kBlue);

    HDSTFiiter DSTFitter(infilelist, false, false, nEvents);
    std::vector<Int_t> pids{ 14, 11, 14, 9 };
    TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);

    DSTFitter.addFitterTask("4c", pids, ppSystem);

    outfile.cd();
    hmLam_prefit->Write();
    hmLam_post4C->Write();
    outfile->Close();
}