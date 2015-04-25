
/*
//	comptonStatus Class
//      keeps track of helicity, mps number, etc. 
//     usses EPICS and auxiliary VME data
//    working with delayed helicity info, etc. should go here
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
#include "comptonParams.h"
#include "bankstructure.h"
#include "fadcdata.h"
#include "vmeauxdata.h"

#ifndef comptonStatus_h
#define comptonStatus_h

#define HELPLUS 1
#define HELMINUS 0

enum laserStateFlag_t {LASER_RIGHT, LASER_LEFT, LASER_RIGHTOFF, LASER_LEFTOFF,
		       LASER_UNKNOWN};


class comptonStatus {
 public:
  comptonStatus();
  comptonStatus(comptonParams* theParamsIn);
  int DefineStatusBranches(TTree *mytree);
  int newRun();  //called at start of run to initizlize status variables
  bool newMPS(int codaEventNumber,fadcdata* theFADCData, vmeauxdata* theVMEauxdata); //called at new MPS
  int GetHelicityState();  //returns -1 if unknown, otherwise 0 or 1 
  int GetCountMPS(){ return countMPS;};
  int GetLaserState(){ return currentLaserState;};
  bool SetLaserState();
  float GetCalibratedBCM(){ return calbcm;};
  bool IsBeamOn(){return beamOn;};
 private:
  bool updateHelicity(int currentHelicityBit); //deals with delayed helicity, etc. Call once per MPS
  bool predictHelicity(int currentHelicityBit); //called by updateHelicity if in delayed mode
  int updateShiftRegister(int hRead);
  comptonParams* theParams;
  int runnum;
  int countMPS;      
  int countPairs;   //count helicity pairs, etc.
  int subcountPairs;  //first or second MPS within helicity pair?
  int countQuads;
  int subcountQuads;  //where are we within a quad (0 through 3) (-1 for not valid)
  int countEpics; //count number of EPICS events encoutnered
  
  int closkscaler;
  int old_clockscaler;

  bool helicityValid;
  int helicityState;
  int helicityDelay;  //filled in at start of run using comptomParams
  int helicityBits;  //# bits in helicity generating shift register
  int helicityStructure; //#of MPS in each set  (4 for "quads")
  uint32_t fgShreg;     //value for helicity sequence algorithm
  int fgNShreg;    //count since fgShreg was reset
  uint32_t helBitHistory;  //keep track of actual bits read
  //current quad patterns for Bit Readout and for Actual Helicity (loest 4 bits0
  //(same thing if not running with helicty bit delayed)
  uint32_t quadPatternCurrentBit;  //current quad bit pattern
  uint32_t quadPatternHelicity;  //actual helcity pattern (predicted pattern)
  uint32_t pickCurrentBit;  //used to select bit in register to compare to current bit

  int countLaserCycles;  //count number complete laser cycles
  int subcountMPSLaserCycles; //count MPS within current laser cycle
  bool laserStateValid;
  laserStateFlag_t laserStateFlag;  //enum labels of laser modes
  int laserStateFound[6];  //tracks which laser states encountered within laser cycle
  laserStateFlag_t previousLaserState;
  laserStateFlag_t currentLaserState;
  laserStateFlag_t knownLaserState;
  laserStateFlag_t lastKnownLaserState;
  int laser_good;
  int laser_on;
  int laser_count;
  int laser_triplet;
  int first_laser_cycle;
  bool good_laser_triplet;
  bool last_good_triplet;
  int first_laser_cycle_arr[3];
  unsigned int laserstate_clock_start;
  float laserstate_clock;  //time into laser state
  int rtcavpow;
  float cavPowerCalibrated;



  // additional status variables for Root Tree
  int mpsCoda;    //mps counter from CODA file
  int mpsSignal;  //state of mpsSignal from TIR
  int mpsScaler;  //mps scaler readout
    int ithr_near;
    int ithr_far;
    float epbcmu3;
    int buflen;
    unsigned int l1a;
    unsigned int mps;

    unsigned int clockscaler;
    int clock_diff;
    int n4before;
    int n4after;
    int n5before;
    int n5after;
    int dithering;
    int dac;
    int rampdelay;
    int inttime;

    float crystalHV;
    float crystalCurrent;
    float horizontalFingerHV;
    float horizontalFingerCurrent;
    float verticalFingerHV;
    float verticalFingerCurrent;
    float tablePosX;
    float tablePosY;
    float coda_deadtime;	// ratio of mps to l1a
    int ihwp_in;
    int wein_right;
    float rtcavpol;

    float ep_s1power;
    float ep_s2power;
    float epbpmX;
    float epbpmY;
    float epbeameng;
    float eptransmit;
    int beam_trip;
    int beam_burp;
    int eppol_burp;
    int HV_trip;
    int rate_fluct;
    int epics_dead;
    int read_helicity;
    int cavpow_burp;
    int trip_mps;
    float future_bcm;
    int HVtrip_mps;
    float future_HV;
    int rate_burp;
    int rate_trip;
    int rate_cut;
    int rate_fluct_mps;
    int rate_trip_mps;
    int epics_mps;
    float future_finger;
    float future_rate;
    int helicity;
    int next_helicity;
    Bool_t helvalid;

    int bad_rtcavpow;
    int rtbcmu3;
    int lastcavpol;
    float epcavpow;
    float eppoldir;
    float eppolpct;

    float lasteppol;

    // test0 thru test 11?
    int test0;
    int test1;
    int test2;
    int test3;
    int test4;
    int test5;
    int test6;
    int test7;
    int test8;
    int test9;
    int test10;
    int test11;


    unsigned int bcmscaler;
    bool beamOn;
    float calbcm;
    /* //UNKNOWN from original fadcanal.h */

    /* unsigned int triggerscaler; */
    /* unsigned int verticalfingerscaler; */
    /* unsigned int horizontalfingerscaler; */
    /* unsigned int old_bcmscaler; */
    /* unsigned int old_clockscaler; */
    /* unsigned int old_triggerscaler; */
    /* unsigned int old_l1a; */
    /* int old_mps;   //allow -1 for initial value */

};


#endif


