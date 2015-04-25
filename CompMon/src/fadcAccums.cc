//*****************************
//*  fadcAccums.cc
//*  created 02/06/2015 G.B. Franklin
//*  For processing Compton FADC accumuator data
//*
//*****************************
#include <iostream>
using namespace std;
#include "bankstructure.h"
#include "fadcAccums.h"


//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fadcAccums::fadcAccums(){
  //constructor
  return;
}
fadcAccums::fadcAccums(comptonParams* theParamsIn){
  //constructor
  theParams=theParamsIn;  //pointer to parameter handeling class
  return;
}
void fadcAccums::newRun(){
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Define histos
int fadcAccums::DefineHistos() {
  float accumFullScale=2E6;
  float scaledRange=accumFullScale/4.0;
  TString label;
  TString title;
  for(int accum=0; accum<8; accum++){
    label=Form("hAcc_acc%d_All",accum);
    title=Form("Accum %d all MPS",accum);
    hAcc_acc_All[accum]=new TH1F(label,title,
			   5000,-accumFullScale,accumFullScale);

    label=Form("hAcc_acc%d_LaserOn",accum);
    title=Form("Accum %d /bcm (Laser On) ",accum);
    hAcc_acc_laserOn[accum]=new TH1F(label,title,
			   5000,-scaledRange,scaledRange);
    label=Form("hAcc_acc%d_LaserOff",accum);
    title=Form("Accum %d /bcm (Laser Off) ",accum);
    hAcc_acc_laserOff[accum]=new TH1F(label,title,
			   5000,-scaledRange,scaledRange);
  }
  //
  return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcAccums::DefineTree(comptonStatus* theComptonStatus){
  // data output to tree for each Compton trigger (summed pulses)
  quadWiseTree=new TTree("quadWise","Accumulator data organized by helicity quadsxs");
  //now add on variables from comptonStatus 
  theComptonStatus->DefineStatusBranches(quadWiseTree);
  return 0;
}
int fadcAccums::DoAccums(vmeauxdata* theVMEauxdata,
			 fadcdata *theFADCdata, comptonStatus* theComptonStatus){
  //
  int laserState=theComptonStatus->GetLaserState();
  int chan=theParams->GetCalorimeterChannel();;  //assume calorimeter on channel 0 of FADC for now
  if(theFADCdata->IsAccumValid(chan)){
    for(int accum=0; accum<6; accum++){
      accraw[accum] = theFADCdata->GetAccumValue(chan, accum); //64 bit
      nacc[accum] = theFADCdata->GetAccumNumSamples(chan,accum);		
    }
    // Now handle accumulator combinations (Accums 6 and 7)
    accraw[6] = accraw[1] + accraw[2];
    accraw[7] = accraw[2] + accraw[3];
    nacc[6] = nacc[1] + nacc[2];
    nacc[7] = nacc[2] + nacc[3];
    int64_t IntegratedPed;
    for(int accum=0; accum<8; accum++){
      IntegratedPed =(theParams->GetPedestal())*nacc[accum];
      accsig[accum]=IntegratedPed -accraw[accum];
    }
    Double_t bcm=theComptonStatus->GetCalibratedBCM();
    Double_t tmp;
    //Fill histograms only if "BeamOn"
    if(theComptonStatus->IsBeamOn()){
      for(int accum=0; accum<8; accum++){
	hAcc_acc_All[accum]->Fill(accsig[accum]);
	tmp=accsig[accum]/bcm;   //Pedistal correct accum divided byh bcm
	if(laserState==LASER_RIGHT || laserState==LASER_LEFT){
	  hAcc_acc_laserOn[accum]->Fill(tmp);
	}else if(laserState==LASER_RIGHTOFF || laserState==LASER_LEFTOFF){
	  hAcc_acc_laserOff[accum]->Fill(tmp);
	}
      }
    }
    quadWiseTree->Fill();
  }
  return 0;
}
  
