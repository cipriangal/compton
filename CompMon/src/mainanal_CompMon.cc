// mainanal_CompMon.cc
//
//  analysis routine for fadccoda
//  Simplified version designed for monitoring Compton fadc data
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <bitset>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "root.h"
#include "TString.h"
#include "THaEpics.h"
#include "THaCodaFile.h"
#include "THaCodaData.h"
#include "THaEpics.h"
#include "TTree.h"
#include "bankstructure.h"         //coda bank structure
//COMPTON fadc classes

#include "comptonStatus.h"         //tracks helicity, laser state, etc.
#include "comptonParams.h"
#include "fadcdata.h"         //unpacks fadc dat
#include "vmeauxdata.h"     //unpacks VME auxillary data
#include "fadcTriggered.h"    //analyzes FADC triggered data
#include "fadcAccums.h"      //analyzes FADC accumulator data

#include "epicsComptonData.h"      //extracts Compton EPICS data
int main(int argc, char** argv)
{
  //+++++++++++++++++++++++++++++++
  //  initialize instances of classes, etc.
  //
  //Hall A general CODA and EPICS classes
  THaCodaFile *codaData;    //the Compton data
  THaEpics* epicsHandler=new THaEpics();  //handles epics data
  //
  // Compton FADC-specific classes
  epicsComptonData* theEPICSdata=new epicsComptonData();//holds Compton epics
  comptonParams* theParams=new comptonParams();   //handels run-dependent  parameters
  fadcdata* theFADCdata=new fadcdata(theParams);  // unpacks FADC data
  //tracks helicty, laser state, etc.
  comptonStatus* theComptonStatus=new comptonStatus(theParams);
  //auxillary data from VME crate
  vmeauxdata* theVMEauxdata=new vmeauxdata(theParams);
//triggered and summed fadc data
  fadcTriggered* theTriggered=new fadcTriggered(theParams);
  //accumulator data
  fadcAccums* theAccums=new fadcAccums(theParams);

  codaData=new THaCodaFile();
  //
  // Setup intput and output files.  
  //
  TString fileName ;
  int run;
  // 
  // Open CODA data file
  if (argc>1)
  {
    run = atoi(argv[1]);
    //    fileName=Form("/home/franklin/HallA/newCompton/data/run%d.dat",run);
    fileName=Form("lnkCompMon.input");  //let script setup link to input file
    if (codaData->codaOpen(fileName, "r") != 0)
      {
	cout << "\nERROR: Cannot open CODA file\n";
	return -1;  
      } 	// end can't open CODA file

    cout<<"input file: "<<fileName<<endl;
  }else{
    cout<< "Error:  Called with no run number arguement"<<endl;
  }
  int maxEvents;
  if(argc>2){
    maxEvents=atoi(argv[2]);
  }else{
    maxEvents=-1;  //flag for no maxEvents count
  }
  printf("DEBUG maxEvents set to %d\n",maxEvents);
  //
  // Open ROOT outputfile  (keep it simple for now)
  //  char* outfilename =Form("fadc.root");
  char* outfilename =Form("lnkCompMon.output");
  cout<<"output file: "<<outfilename<<endl;

  TFile* rootfile = new TFile(outfilename,"RECREATE","fadc data");
  rootfile->SetCompressionLevel(0);
  //
  // setup parameters  (edit fadcparams.h to change values)
  // then setup histograms and trees
  //
  theParams->newRun(run);  //initialize run-dependent parameters
  theComptonStatus->newRun();  //initialize status-tracking 
  theTriggered->DefineHistos();   
  theTriggered->DefineTriggeredTree(theComptonStatus);
  theAccums->DefineHistos();
  theAccums->DefineTree(theComptonStatus);
  // initialize for new run
  int status;
  int eventType;
  uint32_t* pROCdata;    //will get pointer to ROC data for each event
  int counter=0;     //counts all CODA events read
  int accumulatorCounter=0;  //counts physics events (accumulator readouts)
  //  int fadcROCnum=6;   //search for ROC number 6 data 
  int fadcROCnum=0;   //search for ROC number (zero-> use first ROC found) 

  int verbose=0;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // "first pass" analyzes the FADC data, etc.
  //
  //Loop through CODA data
  // codaData->codaClose();
  // codaData->codaOpen(fileName, "r");
  printf("BEGIN PASS ONE\n");
  int calorimeterChan=0;   //probably should get this from params
  while (eventType>=-2 && eventType!=20  &&
	 (counter<maxEvents||maxEvents<0)){
    counter++;
    status=codaData->codaRead();  //read a CODA event buffer in
    if(status<0) break;
    //codaData->evDump(1);
    eventType=codaData->UnpackEventHeader(verbose);
    //1 -> standard ph ysics event
    //16 synch event, 17 PreStart, 18 Go, 19 Pause, 20 End
    // 131 epics event?
    if(eventType<0) {
      printf("event %d eventcode %x", counter, eventType);
    } else if(eventType==1) {
      pROCdata=codaData->getROCdata(fadcROCnum,0); //assume any ROC is THE ROC for now
      if(pROCdata!=NULL){
	theFADCdata->newMPS();     //clear counters, ,etc.
	status=theFADCdata->UnpackAllSubbanks(codaData,theVMEauxdata);
	int codaEventNumber=codaData->GetEventNumber();
	//
	// finished upacking raw data,  now start analysis 
	// determine laserState, helicity, beamON, etc. with theComptonStatus
	//
	status=theComptonStatus->newMPS(codaEventNumber,theFADCdata,theVMEauxdata);

      	if(theFADCdata->IsSumsValid(calorimeterChan)){   //trigged sums
	  status=theTriggered->DoSummedPulses(theVMEauxdata,
					    theFADCdata,theComptonStatus);
	}
	if(theFADCdata->IsSamplesValid(calorimeterChan)){   //waveform samples
	  status=theTriggered->DoSampledWaveforms(theFADCdata,theComptonStatus);
	}
	if(theFADCdata->IsAccumValid(calorimeterChan)){  //FADC accumulator data
	  status=theAccums->DoAccums(theVMEauxdata,theFADCdata,theComptonStatus);
	  accumulatorCounter++;
	}
      } else if(eventType==131){  //EPICS
	theEPICSdata->Unpack(epicsHandler,codaData->getEvBuffer());
	//      theEPICSdata->DebugDump();
      }
      if(counter%1000==0){
	status=theTriggered->DoNormalizedHistos();
      }
      if(counter%10000==0){
	printf(" counter %d eventType %d \n",counter,eventType);
      }
      
      if(eventType==20) break;
    }	// end of first-pass analysis
  }
  rootfile->Write();
  codaData->codaClose();
  rootfile->Close();
  return 0;
}
