/**
 * HDecayBuilder.h
 *
 *
 */

#ifndef HDECAYBUILDER_H
#define HDECAYBUILDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"

using std::cout;
using std::endl;

class HDecayBuilder
{
private:
    std::vector<HRefitCand> fCands;
    int fVerbose;

public:
    HDecayBuilder();
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    void FillData(HParticleCand * cand, HRefitCand & outcand, double arr[], double mass);

};

#endif /* HDECAYBUILDER_H */
