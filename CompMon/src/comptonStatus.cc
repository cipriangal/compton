//*****************************
//*  comptonStatus.cc
//*  created 02/11/2015 G.B. Franklin
//*  Designed for keeping track of Compton status (helicity, laser on etc.)
//*
//*****************************
#include <iostream>
using namespace std;
#include "bankstructure.h"
#include "comptonParams.h"  //defines hardwired but editable parameters
#include "fadcdata.h"
#include "comptonStatus.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
comptonStatus::comptonStatus(){
  //constructor
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
comptonStatus::comptonStatus(comptonParams* theParamsIn){
  //constructor
  theParams=theParamsIn;
  return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::DefineStatusBranches(TTree* mytree){
// External-data branches containing status info shared by all tree
  mytree->Branch("helicityState", &helicityState, "helicityState/I");
  mytree->Branch("laserState", &currentLaserState, "laserstate/I");
  mytree->Branch("beamOn", &beamOn, "beamOn/I");
  mytree->Branch("mpsCoda", &mpsCoda, "mpsCoda/I"); //mps # read from Coda header
  mytree->Branch("mpsSignal", &mpsSignal, "mpsSignal/I"); //mps signal state
  mytree->Branch("mpsScaler", &mpsScaler, "mpsScaler/I");  //scaler mps read
  mytree->Branch("mpsAnalyzed", &countMPS, "mpsAnalyzed/I"); //mps counter
  mytree->Branch("countEpics", &countEpics, "countEpics/I");
  mytree->Branch("countLaserCycles", &countLaserCycles, "countLaserCycles/i");


  mytree->Branch("ithrnear", &ithr_near, "ithrnear/I");
  mytree->Branch("ithrfar", &ithr_far, "ithrfar/I");
  mytree->Branch("epbcmu3", &epbcmu3, "epbcmu3/F");
  mytree->Branch("buflen", &buflen, "buflen/I");
  mytree->Branch("l1a", &l1a, "l1a/i");

  mytree->Branch("clock", &clockscaler, "clock/i");
  mytree->Branch("clockdiff", &clock_diff, "clockdiff/I");
  mytree->Branch("n4before", &n4before, "n4before/I");
  mytree->Branch("n4after", &n4after, "n4after/I");
  mytree->Branch("n5before", &n5before, "n5before/I");
  mytree->Branch("n5after", &n5after, "n5after/I");
  mytree->Branch("dithering", &dithering, "dithering/I");
  mytree->Branch("dac", &dac, "dac/I");
  mytree->Branch("rampdelay", &rampdelay, "rampdelay/I");
  mytree->Branch("inttime", &inttime, "inttime/I");

  mytree->Branch("crystalHV", &crystalHV, "crystalHV/F");
  mytree->Branch("crystal_current", &crystalCurrent, "crystal_current/F");
  mytree->Branch("horizontalFingerHV", &horizontalFingerHV, "horizontalFingerHV/F");
  mytree->Branch("horizontalFingerCurrent", &horizontalFingerCurrent, "horizontalFingerCurrent/F");
  mytree->Branch("verticalFingerHV", &verticalFingerHV, "verticalFingerHV/F");
  mytree->Branch("verticalFingerCurrent", &verticalFingerCurrent, "verticalFingerCurrent/F");
  mytree->Branch("tablePosX", &tablePosX, "tablePosX/F");
  mytree->Branch("tablePosY", &tablePosY, "tablePosY/F");
  mytree->Branch("coda_deadtime", &coda_deadtime, "coda_deadtime/F");
  mytree->Branch("ihwp_in", &ihwp_in, "ihwp_in/I");
  mytree->Branch("wein_right", &wein_right, "wein_right/I");

  //mytree->Branch("rtcavpow", &rtcavpow, "rtcavpow/F");
  mytree->Branch("rtcavpow", &rtcavpow, "rtcavpow/F");
  mytree->Branch("s1power", &ep_s1power, "s1power/F");
  mytree->Branch("s2power", &ep_s2power, "s2power/F");
  //  mytree->Branch("saclayrate", &epsaclay_rate, "saclayrate/F");
  mytree->Branch("bpmx", &epbpmX, "bpmx/F");
  mytree->Branch("bpmy", &epbpmY, "bpmy/F");
  mytree->Branch("bcm", &calbcm, "bcm/F");
  mytree->Branch("beameng", &epbeameng, "beameng/F");
  mytree->Branch("transmit", &eptransmit, "transmit/F");
  mytree->Branch("beam_trip", &beam_trip, "beam_trip/I");
  mytree->Branch("beam_burp", &beam_burp, "beam_burp/I");
  mytree->Branch("eppol_burp", &eppol_burp, "eppol_burp/I");
  mytree->Branch("HV_trip", &HV_trip, "HV_trip/I");
  mytree->Branch("rate_fluct", &rate_cut, "rate_fluct/I");
  mytree->Branch("epics_dead", &epics_dead, "epics_dead/I");
  mytree->Branch("test0", &test0, "test0/I");
  mytree->Branch("test1", &test1, "test1/I");
  mytree->Branch("test2", &test2, "test2/I");
  mytree->Branch("test3", &test3, "test3/I");
  mytree->Branch("test4", &test4, "test4/I");
  mytree->Branch("test5", &test5, "test5/I");
  mytree->Branch("test6", &test6, "test6/I");
  mytree->Branch("test8", &test8, "test8/I");
  mytree->Branch("test9", &test9, "test9/I");
  mytree->Branch("test10", &test10, "test10/I");
  mytree->Branch("test11", &test11, "test11/I");
  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::newRun(){
  //Init all counters and status
  countMPS=0;
  countPairs=0;
  subcountPairs=0;
  countQuads=0;
  subcountQuads=-1;   // 0 through 3 for valid MPS within a Quad
  countLaserCycles=0;
  subcountMPSLaserCycles=0; //number MPS within current laser mode
  helicityValid=false;
  helicityState=2;
  fgShreg=0;   //helcity register decoding
  fgNShreg=0;   //counts since fgShreg reset
  helBitHistory=0;  //shift register for actual helicity bit read
  laserStateValid=false;
  currentLaserState=LASER_UNKNOWN;  //enum
  for(int i=0; i<6; i++){
    laserStateFound[i]=false;
  }
  //look up run-dependent helicity info
  helicityDelay=theParams->GetHelicityDelay(); //run dependent
  helicityBits=theParams->GetHelicityBits();
  helicityStructure=theParams->GetHelicityStructure();
  printf("Initializing Run      %d\n",theParams->GetRun());
  printf("   Helicity Delay     %d\n", helicityDelay);
  printf("   Helicity Bits      %d\n", helicityBits);
  printf("   Helicity Structure %d\n",helicityStructure);
  return 0;
}
bool comptonStatus::newMPS(int codaEventNumber, fadcdata* theFADCdata, vmeauxdata* theAuxData){
  countMPS++;   //mps count via count of analyzed Event 1s
  mpsCoda=codaEventNumber;  //mps count via CODA Event 1 header
  mpsScaler=theAuxData->GetMPSScaler();  //mps count via VME scaler
  mpsSignal=theAuxData->GetMPSSignal();
  dithering=theAuxData->GetDithering();
  rtcavpol= theAuxData->GetCavityPolBit();  
  int cavityPowerBit=theAuxData->GetCavityPowerBit(); 
  helicityValid=updateHelicity(theAuxData->GetHelicityBit());  //deals with delayed helicity, etc.

  //laser cavity pwoer
  old_clockscaler=clockscaler;
  clockscaler = theAuxData->GetClockScaler();
  clock_diff=clockscaler-old_clockscaler;
  if(clock_diff<0){
    printf("\n WARNING old clock= %d, new  clock=%d",old_clockscaler,clockscaler);
  }else{
    cavPowerCalibrated=( 4.0E7*theParams->GetCavPowCalibration()*
    theAuxData->GetCavityPowerScaler() )/clock_diff;
  }
  bcmscaler =theAuxData->GetGatedBCM();  //bcm value from VTF gated scaler
  calbcm = (float) bcmscaler/ (float)clock_diff;
  calbcm*= 4E7*theParams->GetBCMCalibration();
  beamOn=calbcm > (theParams->GetHighBCMCutoff());
  /*
  printf("mpsScaler %d\n",mpsScaler);
  printf("mpsSignal %d\n",mpsSignal);
  printf("dithering %d\n",dithering);
  printf("rtcavpol  %f\n",rtcavpol);
  printf("cavityPowerBit %d\n",cavityPowerBit);
  printf("rtcavPower  %f\n",rtcavpow);
  printf("cavityPowerScaler %d\n",theAuxData->GetCavityPowerScaler());
  printf("clock_diff %d\n",clock_diff);
  printf("CavCalibration %f\n",theParams->GetCavPowCalibration());
  // determine laser state, where we are in a laser triplet, etc.
  */
  bool transition=SetLaserState(); //determine laser state
  //transfer in scalers.
  test0 = theAuxData->GetIPScaler(0);
  test1 = theAuxData->GetIPScaler(1);
  test2 = theAuxData->GetIPScaler(2);
  test3 = theAuxData->GetIPScaler(3);
  test4 = theAuxData->GetIPScaler(4);
  test5 = theAuxData->GetIPScaler(5);
  test6 = theAuxData->GetIPScaler(6);
  test7 = theAuxData->GetIPScaler(7);
  test8 = theAuxData->GetIPScaler(8);
  test9 = theAuxData->GetIPScaler(9);
  test10= theAuxData->GetIPScaler(10);
  test11= theAuxData->GetIPScaler(11);

  return transition;
}
int comptonStatus::GetHelicityState(){
  if(helicityValid) {
    return helicityState;
  }else{
    return -1;
  }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//e
bool comptonStatus::updateHelicity(int currentBit){
  // updates Helicity and quad status
  //
  //
  //QuadA bits 0110 = 0x06
  //QuadB bits 1001 = 0x09
  helBitHistory= ( helBitHistory<<1 | currentBit) & 0xFF;
  //update subQuad and Quad count if we've been synched
  if(subcountQuads>=0){
    subcountQuads++;
    if(subcountQuads>3) {
      subcountQuads=0;
      countQuads++;
    }
  }
  //if not synched to quads, look at last 8 MPSs to see if a transition
  if(subcountQuads<0){
    if(helBitHistory==0x69 || helBitHistory==0x96) {
      //QuadA to Quad B transition or QuadB to QUadA transition
      subcountQuads=3; //this was last MPS within a quad
      printf("synched quads\n");
    }
  }
  if(helicityDelay==0){
    helicityState=currentBit;  //nothing fancy if helicity bit is not delayed
    helicityValid=true;
  }else if (helicityDelay==8&&helicityBits==30&&helicityStructure==4) {
      helicityValid=predictHelicity(currentBit);
  }else{
    printf("WARNING: unsported helcity structure\n");
    helicityValid=false;
  }
  return helicityValid;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
bool comptonStatus::predictHelicity(int currentBit){
    // uses first events to  train shift register, then predicts quad pattern
    // in the future
  //This version hardwires to 30 bit shift register, 8 MPS delay (q qauds)
  const int pickBit[4]={1,2,4,8};  //used to pick bit within subquad
  int bitCheck;
  int futureBit;
  bool predictValid;
  //
  //first, update shiftregister if at the start of a quad
  // if subcountQuad<30, we are still training the shift regiseter
  //  when subcountQuads==30, we update shiftregister to get prediction'
  // of future quad pairs
  //  fgShreg used by "updateShiftRegister" to generatre psuedo-randoms
  //  fgNShreg used to count training bits (need 30 to fill register
  //   pickCurrentBit will have 1 bit set to pick off fgShreg bit corresponding
  //     to the prediction of the current Bit
  //  quadPatterCurrentBit  holds the 4-bit patter of the predicted current quad
  //
  if(subcountQuads==0){
    if(fgNShreg<30){     //are we still training the shift register?
      futureBit = updateShiftRegister(currentBit); //yes, feed it the bit
      fgNShreg++;
      bitCheck=-1;  //nothing to check yet
      if(fgNShreg==30) {
	//finished training.  jump forward in in predictions and
	//also record bit within shiftregister for "currentBit" preddiciton"
	for(int i=0; i<helicityDelay/4; i++){ //finished training, jump ahead 
	  futureBit = updateShiftRegister(2);
	}
	pickCurrentBit=1;
	pickCurrentBit=pickCurrentBit<<helicityDelay/4; 
      } else {
	futureBit=-1;  //not able to predict the future yet
	pickCurrentBit=1;
      }
    }else{
      futureBit = updateShiftRegister(2); //no, let it decide
      bitCheck=(fgShreg &pickCurrentBit)>0; //check in register for prediction of currrent bit
    }
    if( bitCheck >0){   //check 8 bits back in register prediction of current bit
      quadPatternCurrentBit=0x09;
    }else{
      quadPatternCurrentBit=0x06;
    }
    if( futureBit >0) { //check predicted quad
      quadPatternHelicity=0x09;
    }else{
      quadPatternHelicity=0x06;
    }
  }
//see currentBit agrees with predicted currentbit
  if(fgNShreg>=30 && subcountQuads>=0){
    bitCheck= (quadPatternCurrentBit&pickBit[subcountQuads])>0;
    if(bitCheck!=currentBit){
      printf("helicity bit check failed\n");
      fgNShreg=0;   //assume shift register needs retaining
      subcountQuads=-1;  //assume we need to resynch quads
    }
  }
  if(fgNShreg>=30 && subcountQuads>=0){
    helicityState=(quadPatternHelicity&pickBit[subcountQuads])>0;
    predictValid=true;
  }else{
    helicityState=2;
    predictValid=false;
  }
  return predictValid;
} 
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
int comptonStatus::updateShiftRegister(int hRead){
  // updates helcity-tracking shifdt-register used for pseudo-randoms
  //  input hRead= current helicty bit read for training registrer
  //  input hRead= 2 to return prediction of next bit
  //30 bit structure
  uint32_t bit7    = (fgShreg & 0x00000040) != 0;
  uint32_t bit28   = (fgShreg & 0x08000000) != 0;
  uint32_t bit29   = (fgShreg & 0x10000000) != 0;
  uint32_t bit30   = (fgShreg & 0x20000000) != 0;
  
  int newbit = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;
  fgShreg = ( (hRead == 2 ? newbit : hRead) | (fgShreg << 1 )) & 0x3FFFFFFF;
  return newbit;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Use TIR and EPICS data to figure out current/previous laser state
  bool comptonStatus::SetLaserState(){
  bool transitioned = false;
  //get threshold parameters used for cavity power
  int rtcavpowon = theParams->GetCavPowOn();
  int rtcavpowoff = theParams->GetCavPowOff();

  // Set previous laser state
  previousLaserState = currentLaserState;
  lastKnownLaserState = knownLaserState;
  
  // Set new current laser state
  // Due to the unreliability of rtcavpow, we are eliminating it from the test
  if(cavpow_burp){
    currentLaserState=LASER_UNKNOWN;
  }
  else if (cavPowerCalibrated < rtcavpowoff){   // Cavity is off!
      if(rtcavpol > 0.5){
      currentLaserState = LASER_LEFTOFF;
    } else {
      currentLaserState = LASER_RIGHTOFF;
    }
  } else if (cavPowerCalibrated < rtcavpowon ){ 
    currentLaserState = LASER_UNKNOWN; // intermediate cavity power
  } else {
    if ( rtcavpol > 0.5){ 
      currentLaserState = LASER_LEFT;
    } else {
      currentLaserState = LASER_RIGHT;
    }
  }										    
  if(currentLaserState == LASER_LEFT || currentLaserState == LASER_RIGHT){
    laser_on=1;
  }
  else if(currentLaserState == LASER_LEFTOFF ||
	  currentLaserState == LASER_RIGHTOFF){
    laser_on=0;
  }
  else if(currentLaserState == LASER_UNKNOWN){
    laser_on=-1;
  }

  // Now let's compare to the previous laser state

  if(currentLaserState != LASER_UNKNOWN){
  knownLaserState = currentLaserState;
  }
  if (knownLaserState != lastKnownLaserState ){
    transitioned = true;
    laser_good = 0;
    laserstate_clock_start=clockscaler;
    laser_count++;
    laser_triplet++;
  } else {
    laserstate_clock= (float)(clockscaler-laserstate_clock_start)/4.e7;
  }
  if(transitioned){
    if(laser_triplet==3||(first_laser_cycle&&!laser_on&&laser_triplet!=2)){
      laser_triplet=0;
    }
    //cerr<<"laser triplet "<<laser_triplet<<" laserstate "<<current_laser_state<<endl;
    first_laser_cycle_arr[laser_triplet]=laser_on;
    last_good_triplet = good_laser_triplet;
    good_laser_triplet=false;
    if(laser_triplet==2){
      if(first_laser_cycle_arr[0]==0&&first_laser_cycle_arr[1]==1&&first_laser_cycle_arr[2]==0){
        good_laser_triplet=true;
      }

      if((!good_laser_triplet||lastcavpol==rtcavpol)&&!first_laser_cycle){
        //if((!good_laser_triplet)&&!first_laser_cycle){
	
	// InitLaserCounters();
	
        first_laser_cycle=1;
        cerr<<"Error, problem with laser_triplet id"<<endl;
        cerr<<"laser on "<<first_laser_cycle_arr[0]<<" "<<first_laser_cycle_arr[1]<<" "<<first_laser_cycle_arr[2]<<" laser triplet "<<laser_triplet<<" laserstate "<<currentLaserState<<endl;
      }
      for(int i=0; i<3; i++){
        first_laser_cycle_arr[i] = -1;
      }
      }
      if(good_laser_triplet&&lastcavpol!=rtcavpol&&first_laser_cycle){
        first_laser_cycle=0;
        //cerr<<"first laser cycle done"<<endl;
      }

  }
  return transitioned;
  }

