//*****************************
//*  fadcTriggered.cc
//*  created 02/06/2015 G.B. Franklin
//*  Designed for processing Compton FADC triggered data
//*
//*****************************
#include <iostream>
using namespace std;
#include "bankstructure.h"
#include "fadcTriggered.h"

#define TIMEWINDOWlow 566
#define TIMEWINDOWhi 150

#define TIMEWINDOW2low 150
#define TIMEWINDOW2hi 500

#define SUMWINDOWlow 10000
#define SUMWINDOWhi 30000

#define SNAPLIMIT 1     //-&-

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fadcTriggered::fadcTriggered(){
  //constructor
  return;
}
fadcTriggered::fadcTriggered(comptonParams* theParamsIn){
  //constructor
  theParams=theParamsIn;  //pointer to parameter handeling class
  return;
}
void fadcTriggered::newRun(){
  mpsLaserOnCount=0;
  mpsLaserOffCount=0;
  bcmLaserOnSum=0.;
  bcmLaserOffSum=0.;
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Define histos
int fadcTriggered::DefineHistos() {
  int fullscale=100000;
  int offset=-1000;
  //  int smallscale=5000;
  //  int smalloffset=-1000;
  hTrig_numSums=new TH1F("hTrig_numSums","Number Summed FADC Triggers in MPS",
			   100,0,100);
  hTrig_numSamples=new TH1F("hTrig_numSamples",
			    "Number of WF samples in MPS",100,0,100);
  //
  // histograms of integrated triggered pulses
  hTrig_sums_All=new TH1F("hTrig_sums_All","Sums of All Triggered Pulses",
			 5000,offset,fullscale+offset);

  hTrig_sums_laserOn=new TH1F("hTrig_sums_laserOn",
			      "Sums of Laser On Triggered Pulses: Laser On",
			      5000,offset,fullscale+offset);
  hTrig_sums_laserOff=new TH1F("hTrig_sums_laserOff",
			      "Sums of Laser On Triggered Pulses: Laser Off",
			      5000,offset,fullscale+offset);
  hTrig_sums[0]=new TH1F("hTrig_sums_L_P","Sums for (Laser Left)-(Helicty +1)",
			 5000,offset,fullscale+offset);

  hTrig_sums[1]=new TH1F("hTrig_sums_L_N","Sums for (Laser Left)-(Helicty -1)",
			 5000,offset,fullscale+offset);
  hTrig_sums[2]=new TH1F("hTrig_sums_R_P","Sums for (Laser Right)-(Helicty +1)",
			 5000,offset,fullscale+offset);
  hTrig_sums[3]=new TH1F("hTrig_sums_R_N","Sums for (Laser Right)-(Helicty -1)",
			 5000,offset,fullscale+offset);
  //
  // Histograms that count triggers sorted by Laser State
  hTrig_MPS=new TH1F("hTrig_MPS",
			       "Num MPS  vs Laser State",10,0,10);
  labelLaserSortedHistos(hTrig_MPS);


  //
  hTrig_MPS_BeamOn=new TH1F("hTrig_MPS_BeamOn",
			       "Num MPS with Beam On  vs Laser State",10,0,10);
  labelLaserSortedHistos(hTrig_MPS_BeamOn); //set labels of x-axis
 //
  hTrig_BCM=new TH1F("hTrig_BCM", "BCM  vs Laser State",10,0,10);
  labelLaserSortedHistos(hTrig_BCM);  //set labels of x-axis
  //
  hTrig_Trigs_Accepted=new TH1F("hTrig_Trig_.Accepted",
				 "Num Trig Evnts Accepted vs Laser State",10,0,10);
 labelLaserSortedHistos(hTrig_Trigs_Accepted);
  //
  hTrig_Trigs_Scaler=new TH1F("hTrig_Trigs_Scaler",
				 "Num Triggers vs Laser State",10,0,10);
 labelLaserSortedHistos(hTrig_Trigs_Scaler);
  //
  //normalized histos
  hNorm_Trigs_Scaler=new TH1F("hNorm_Trigs_Scaler",
				 "Num Triggers per MPS vs Laser State",10,0,10);
  labelLaserSortedHistos(hNorm_Trigs_Scaler);

  hNorm_sums_subtracted=new TH1F("hNorm_sums_subtracted",
			      "Background Subtracted Triggered Pulses",
			      5000,offset,fullscale+offset);
  hNorm_sums_asym=new TH1F("hNorm_sums_asym",
			      "Energy-Dependent Asymmetry",
			      5000,offset,fullscale+offset);
  //
  //waveform snapshop histos
  hTrig_wf=new TH1F("hTrig_wf","Sum of Sampled Waveforms",500,0,500);
  return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DefineTriggeredTree(comptonStatus* theComptonStatus){
  // data output to tree for each Compton trigger (summed pulses)
  triggerWiseTree=new TTree("triggerwise","Pulse-wise triggered data");
  triggerWiseTree->Branch("sum",&sumVal,0);
  //now add on variables from comptonStatus 
  theComptonStatus->DefineStatusBranches(triggerWiseTree);
  //
  // data for sampled waveforms
  snapshotsTree=new TTree("snapshots","sampled snapshops");
  snapshotsTree->Branch("numSamples",&NumSamples,"numSamples/I");
  snapshotsTree->Branch("snap",&snapshot,"snapshot[numSamples]/F");
  //now add on variables from comptonStatus 
  theComptonStatus->DefineStatusBranches(snapshotsTree);
  return 0;
}
void fadcTriggered::labelLaserSortedHistos(TH1F* histo){
  //helper function to DefineHisos
  //
  TAxis* pAxis= histo->GetXaxis();
  pAxis->SetBinLabel(1,"LN");
  pAxis->SetBinLabel(2,"LP");
  pAxis->SetBinLabel(3,"RN");
  pAxis->SetBinLabel(4,"RP");
  pAxis->SetBinLabel(5,"Roff");
  pAxis->SetBinLabel(6,"Loff");
  pAxis->SetBinLabel(7,"Unk");
  return;
}

int fadcTriggered::DoSummedPulses(vmeauxdata* theVMEauxdata,
				  fadcdata *theFADCdata,
				  comptonStatus* theComptonStatus){
  // histogram triggered data (pre-summed by CODA )
  int chan=theParams->GetCalorimeterChannel(); //fadc channel used for calor.
  if(theFADCdata->IsSumsValid(chan)){
    int originalPedCorr=theFADCdata->GetSumsPedestalSubtracted(chan);
    int numInSum=theFADCdata->GetNumberSamplesSummed(chan);
    double PedCorrection=numInSum*theParams->GetPedestal() - originalPedCorr;
    int numTriggersAccepted=theFADCdata->GetSumsNumberTriggersSummed(chan);
    int numTriggers=theVMEauxdata->GetTriggerScaler();

    mpsCount=theComptonStatus->GetCountMPS();
    helicityState=theComptonStatus->GetHelicityState();
    laserState=theComptonStatus->GetLaserState();
    float bcm=theComptonStatus->GetCalibratedBCM();
    bool beamOn= theComptonStatus->IsBeamOn();
    int sort;  
    if(laserState==LASER_LEFT){
      if(helicityState==0){
	sort=0;
      }else{
	sort=1;
      }
    }else if(laserState==LASER_RIGHT){
      if(helicityState==0){
	sort=2;
      }else{
	sort=3;
      }
    }else if (laserState==LASER_LEFTOFF){
      sort=4;
    }else if (laserState==LASER_RIGHTOFF){	
      sort=5;
    }else{
      sort=6;
    }
    //prepare to ignore LASER_UNKNOWN data
    bool laserOn= (laserState==LASER_RIGHT||laserState==LASER_LEFT);
    bool laserOff= (laserState==LASER_RIGHTOFF||laserState==LASER_LEFTOFF);

    hTrig_MPS->Fill(sort);    //MPS counts by laser cycle
    //following sorted fills are only for BEAM ON condiation
    if(beamOn){
      if(laserOn){
	mpsLaserOnCount++;
	bcmLaserOnSum+=bcm;
      }else if (laserOff){
	mpsLaserOffCount++;
	bcmLaserOffSum+=bcm;
      }
      hTrig_MPS_BeamOn->Fill(sort);
      hTrig_Trigs_Scaler->Fill( sort,numTriggers);
      hTrig_Trigs_Accepted->Fill(sort,numTriggersAccepted);
      hTrig_numSums->Fill(numTriggersAccepted);
      hTrig_BCM->Fill(sort,bcm); //
      for(int i=0; i<numTriggersAccepted; i++){
	sumVal=theFADCdata->GetSums(chan,i)+PedCorrection;  //move to roottree slot
	hTrig_sums_All->Fill(sumVal);
	if(laserOn) hTrig_sums_laserOn->Fill(sumVal);
	if(laserOff) hTrig_sums_laserOff->Fill(sumVal);
	if(sort>=0 && sort<4)  hTrig_sums[sort]->Fill(sumVal);
      }
    }
    triggerWiseTree->Fill();
  }else{
    printf("Invalid Summed Triggered Data\n");
  }
  return 0; 
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DoNormalizedHistos(){
  // calculates normalized histos (from triggered data)
  //really stupid way to copy histogram.  (Must be a better way)
  hNorm_Trigs_Scaler->Add(hTrig_Trigs_Scaler,hTrig_Trigs_Scaler,1.,0.);
  hNorm_Trigs_Scaler->Divide(hTrig_MPS); //normalizer per MPS for each laser state

  //
  if(bcmLaserOffSum>0.0){
    Double_t C1=1.0;
    Double_t C2= -bcmLaserOnSum/(bcmLaserOffSum);  //note negative sign
    hNorm_sums_subtracted->Add(hTrig_sums_laserOn,hTrig_sums_laserOff,C1,C2);
    // build asymmetry histogram
    //start by building denominator
    hNorm_sums_asym->Add(hTrig_sums[1],hTrig_sums[2],1.,1.);
    hNorm_sums_asym->Add(hTrig_sums[0],-1.);  //note negative signs
    hNorm_sums_asym->Add(hTrig_sums[3],-1.);
    //divide by background subtracted sums
    hNorm_sums_asym->Divide(hNorm_sums_subtracted);
  }
  return 0;
}

//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DoSampledWaveforms(fadcdata *theFADCdata, 
				      comptonStatus* theComptonStatus){
  // histogram triggered data (presummed by CODA )
   if(theFADCdata->IsSamplesValid(0)){
    int chan=0;   //only FADC channel zero for now
 //number FADC channel read out
    int TotalSamples=theFADCdata->GetNumberSamples(chan);
 //# samples stored for each sample trigger
    NumSamples=theFADCdata->GetSamplesPerEvent(chan);
    //# of waveforms stored
    int NumEvents= theFADCdata->GetNumberEvents(chan);
    if(TotalSamples!= NumSamples*NumEvents){
      NumEvents=TotalSamples/NumSamples;   //fix if data overflow
      printf("Total Samples doesn't agree with number of events\n");
    }
    mpsCount=theComptonStatus->GetCountMPS();
    helicityState=theComptonStatus->GetHelicityState();
    laserState=theComptonStatus->GetLaserState();
    hTrig_numSamples->Fill(NumEvents);
    pulsestructure Pulse;
    int* bits;
    int* data;
    for(int event=0; event<NumEvents; event++){
      Pulse=theFADCdata->GetPulse(chan,event);
      bits=Pulse.UserBits;
      data=Pulse.Data;
      for(int i=0; i<NumSamples; i++){
 	snapshot[i]=data[i];
	hTrig_wf->Fill(i,data[i]);  //use weighted fill
      }
      snapshotsTree->Fill();
    }
  }else{
    printf("Invalid Waveform Data\n");
  }
  return 0;   
}

