#include <stdio.h>
#include <TROOT.h>
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "bankstructure.h"
#include "comptonParams.h"

#ifndef vmeaux_h
#define vmeaux_h

// Add scaler data?
class vmeauxdata {
  private:
  comptonParams* theParams; // pointer to parameter handling class
  // header data
  int numAuxWords;
    // TIR data
  int helBit;   //helicity bit
  int dithering; //dithering on bit
    int qrt;
    int cavpowBit;  //cavity power Bit (but also look at real-time cav power
    int cavpolBit;
    int evlen;
    int mpsSignal;

    // HAPPEX Timing Board data
    int rampDelay;
    int intTime;
    int DACsetting;   //DAC setting for pulser
    //
    unsigned int scalers[16];
    unsigned int IPscalers[16];
    //private member functions
    void init();  //used to clear private variables
  public:

    vmeauxdata();		        	// Constructor
    vmeauxdata(comptonParams* theParams);  // Constructor

    unsigned int GetMPSScaler()	{return scalers[5];};
    unsigned int GetL1A() 	{return scalers[4];};
    int GetMPSSignal()		{return mpsSignal;};
    int GetHelicityBit()	{return helBit;};

    int GetDithering()  {return dithering;};
    int GetQRT()		{return qrt;};
    int GetCavityPowerBit()	{return cavpowBit;};
    int GetCavityPolBit()	{return cavpolBit;};
    int GetRampDelay()  	{return rampDelay;};
    int GetIntTime()    	{return intTime;};
    int GetDACSetting()         {return DACsetting;};
    unsigned int GetScaler(int i){return scalers[i];};
    unsigned int GetBCMScaler()  {return scalers[0];};
    unsigned int GetClockScaler(){return scalers[3];};
    //    unsigned int GetTriggerScaler(){return scalers[6];};
    unsigned int GetTriggerScaler(){return IPscalers[12];};

    unsigned int GetIPScaler(int i){return IPscalers[i];};
    unsigned int GetGatedBCM(){return IPscalers[0];};
    unsigned int GetGatedClock(){return IPscalers[11];};
    unsigned int GetVtoFBPM2AypIPScaler(){return IPscalers[1];};
    unsigned int GetVtoFBPM2AxmIPScaler(){return IPscalers[2];};
    unsigned int GetVtoFBPM2AxpIPScaler(){return IPscalers[3];};
    unsigned int GetVtoFPowLeftIPScaler(){return IPscalers[5];};
    unsigned int GetVtoFPowRightIPScaler(){return IPscalers[6];};
    unsigned int GetVtoFPowIPScaler(){return IPscalers[7];};
    unsigned int GetVtoFBPM2BymIPScaler(){return IPscalers[8];};
    unsigned int GetVtoFBPM2BypIPScaler(){return IPscalers[9];};
    unsigned int GetVtoFBPM2BxmIPScaler(){return IPscalers[10];};
    unsigned int GetVtoFBPM2BxpIPScaler(){return IPscalers[11];};

    unsigned int GetTriggerIPScaler(){return IPscalers[12];};
    unsigned int GetVerticalFingerIPScaler(){return IPscalers[14];};
    unsigned int GetHorizontalFingerIPScaler(){return IPscalers[15];};


    unsigned int GetTest5(){return IPscalers[5];};
    unsigned int GetTest6(){return IPscalers[6];};
    unsigned int GetCavityPowerScaler(){return IPscalers[7];};

    int Unpack(bankstructure bank);

    // Which module information to dump?	0 = all
    //										Multiples of 2 = Trigger Interface
    //										Multiples of 3 = Timing Board
    //										Multiples of 5 = Scalers
    // e.g. 6 gives you both trigger interface and timing board
    int DebugDump(int whichmodule=0);
};

#endif 	//vmeaux_h
