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

Int_t analysis_user(TString infileList="/lustre/hades/user/rlalik/hades/pp45/sim/simdst/out/pluto_chan_060_events_50000_seed_00*_1_dst_p4500p.root", Int_t nEvents=-1){

    HDSTFiiter DSTFitter(infilelist, false, false, nEvents);
    std::vector<Int_t> pids{ 14, 11, 14, 9 };
    TLorentzVector lv = ();

    DSTFitter.addFitterTask("4c", pids, lv);
}