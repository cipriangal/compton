
/*
//	fadcTriggered
//      Analysis of triggered fadc dada
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

#ifndef fadcTriggeredSums_h
#define fadcTriggeredSums_h

#define MAX_DAC_VALUES 20
#define HELPLUS 1
#define HELMINUS 0

// laserstate now defined in class comptonStatus
//enum laserstate {RIGHT, LEFT, OFF, UNKNOWN};
//enum laserstate {RIGHT, LEFT, RIGHTOFF, LEFTOFF, UNKNOWN};

class fadcTriggered {
 public:
  fadcTriggered();
  fadcTriggered(comptonParams* theParams);
  void newRun();   //init counters at start of a run
  int DefineHistos();
  int DefineTriggeredTree(comptonStatus* theComptonStatus);
  int DoSummedPulses(vmeauxdata* theVMEauxdata,
		     fadcdata *theFADCdata, comptonStatus* theComptonStatus); 
  int DoNormalizedHistos();
  int DoSampledWaveforms(fadcdata *theFADCdata,comptonStatus* theComptonStatus); 
 private:
  comptonParams* theParams;   //pointer to parameter 
  void labelLaserSortedHistos(TH1F* histo); //used to label x-axis of laser-sorted histos
  int NumSamples;  //Filled with number of valid waveform samples (snapshots)
  //
  //Variables output to triggeredWise root tree
  TTree* triggerWiseTree;
  double sumVal;   //fadc triggered data summed over single pulse
  //
  // special tree for snapshots  (sampled waveforms)
  TTree* snapshotsTree;
  int mpsCount;
  int mpsLaserOnCount;
  int mpsLaserOffCount;
  Float_t bcmLaserOnSum;
  Float_t bcmLaserOffSum;
  int helicityState;
  int laserState;
  Float_t snapshot[1000]; //will use this for a snapshot array
  //
  // Triggered Sums data histograms
  TH1F* hTrig_numSums;  //#triggers actually summed each MPS
  TH1F* hTrig_sums_All;  //histo of all summed pulses
  TH1F* hTrig_sums_laserOn;  //histo of laser-on summed pulses
  TH1F* hTrig_sums_laserOff;  //histo of laser-on summed pulses
  TH1F* hTrig_sums[4];  // pulse sums for Laser Left, Positive Helicity Beam. etc


  TH1F* hTrig_MPS;      //count MPSs sorted by Laser State
  TH1F* hTrig_MPS_BeamOn;      //count MPSs with BeamOn sorted by Laser State
  TH1F* hTrig_BCM;    
//BCM sorted by Laser State
  TH1F* hTrig_Trigs_Accepted;  //count accepted triggers sorted by Laser State
  TH1F* hTrig_Trigs_Scaler;   //count trigger (via scaler) sorted by Laser State
  // calculated periodically from other histograms 
 //Normalized histos
  TH1F* hNorm_Trigs_Scaler;   //triggers (via scaler) per MPS
  TH1F* hNorm_sums_subtracted; //background subtracted (laserOn-laserOff)
  TH1F* hNorm_sums_asym;   //asymmetry vs energy
  //Triggered Sampled waveform histograms
  TH1F* hTrig_numSamples;
  TH1F* hTrig_wf;
};

#endif
