
/*
//	fadcAccums
//      Analysis of fadc accumulator dada
 */

#include <stdlib.h>
#include <TROOT.h>
#include "TTree.h"
#include "TMath.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "THaEpics.h"
#include "THaCodaFile.h"
#include "bankstructure.h"
//Compton fadc-specific classes
#include "comptonStatus.h"
#include "comptonParams.h"
#include "fadcdata.h"
#include "vmeauxdata.h"

#ifndef fadcAccums_h
#define fadcAccums_h

#define MAX_DAC_VALUES 20
#define HELPLUS 1
#define HELMINUS 0

// laserstate now defined in class comptonStatus
//enum laserstate {RIGHT, LEFT, RIGHTOFF, LEFTOFF, UNKNOWN};

class fadcAccums {
 public:
  fadcAccums();
  fadcAccums(comptonParams* theParams);
  void newRun();   //init counters at start of a run
  int DefineHistos();
  int DefineTree(comptonStatus* theComptonStatus);
  int DoAccums(vmeauxdata* theVMEauxdata,
		     fadcdata *theFADCdata, comptonStatus* theComptonStatus); 

  int DoNormalizedHistos();

 private:
  comptonParams* theParams;   //pointer to parameter 
  //
  //Variables output to triggeredWise root tree
  TTree* quadWiseTree;
  //accumulator data
  int nacc[8];
  int64_t accraw[8];  //raw acumuilator data
  double accsig[8];   //pedestal corrected accumulator data
  //
  // Accumulator data histograms
  TH1F* hAcc_acc_All[8];  //#triggers actually summed each MPS
  TH1F* hAcc_acc_laserOn[8];  //histo of laser-on summed pulses
  TH1F* hAcc_acc_laserOff[8];  //histo of laser-off summed pulses
};

#endif
