#ifndef HCOVARIANCEKINFIT_H
#define HCOVARIANCEKINFIT_H

//#include "hparcond.h"    is this useful here?

//class HCovarianceKinFit : public HParCond {
class HCovarianceKinFit {
protected:
	TString fSetup;
	bool fMomDepErrors = false;
	Int_t fVerbose = 0;
/*
// for protons
  Int_t    nThetaReg;            // number of theta regions 
  Double_t thetaStep;            // size of one theta region, degrees
  Double_t thetaMiddle0;         // middle of first theta region, degrees
  Double_t momInt1;              // end of first momentum interval, MeV
  Double_t momInt2;              // end of second momentum interval, MeV
  Int_t    nParams;              // number of parameters for one theta region
  Double_t parMomCorrH[264];     // parameters. Size=nThetaReg*nParams
  Short_t  typePar;              // type of parametrisation 

//for electrons  
  Int_t    nParamsElect;         // number of parameters
  Double_t parMomCorrElect[13];   // parameters
  */
  static   HCovarianceKinFit* gCovariance;
  
public:
  HCovarianceKinFit(const Char_t* name    = "CovarianceKinFit",
                     const Char_t* title   = "Covariance matrix estimation for kin fit",
                     const Char_t* context = "CovarianceMatrixEstimate");
  ~HCovarianceKinFit(void) {}
  static HCovarianceKinFit* getObject(void) {return gCovariance;}
  
  void estimateCov(Int_t pid, Double_t mom, double (&covariance)[5]);
  void setSetup(TString run){ fSetup=run; }
  void setMomDepErrors(bool val){ fMomDepErrors=val; }
  /*
  void     clear(void);
  void     putParams(HParamList*);
  Bool_t   getParams(HParamList*);

  Bool_t   setDefaultPar(TString run);
  Double_t getDeltaMom(Int_t pId, Double_t mom, Double_t theta) const; // units: MeV for mom, degrees for theta
  Double_t getCorrMom(Int_t pId, Double_t mom, Double_t theta) const {return mom+getDeltaMom(pId,mom,theta);}
*/
private:
/*
  void     fillParMomCorrH(Int_t size,Double_t *par);
  Double_t getDeltaMomT1(Int_t pId, Double_t mom, Double_t theta) const;
  Double_t getDeltaMomT12(Int_t pId, Double_t mom, Double_t theta) const;
  Double_t getDeltaMomT1419(Int_t pId, Double_t mom, Double_t theta) const;
  Double_t binInter(Double_t b,Double_t rs, const Double_t *par) const;
  */  
  ClassDef(HCovarianceKinFit,1) // Parameter container for energy loss correction. Wozu?
};

#endif  /*!HCOVARIANCEKINFIT_H */
