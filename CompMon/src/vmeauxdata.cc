//*******************************
//    vmeauxdata.cc
//    created 11/15/2008 D. Parno
//    handles VME data from auxiliary modules:
//		HAPPEX timing board
//		Trigger Interface Register
//		16-channel Caen v560 scaler
//*******************************

#include <iostream>
using namespace std;
#include "vmeauxdata.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
vmeauxdata::vmeauxdata(){               // Constructor
  init();  //clear private variables
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
vmeauxdata::vmeauxdata(comptonParams* theParamsIn){   // Constructor
  theParams=theParamsIn;
  init();   //initialize private variables
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void vmeauxdata::init(){
  helBit = -1;
  mpsSignal = 0;
  qrt = 0;
  cavpowBit = 0;
  cavpolBit = 0;
  evlen = 0;
  rampDelay = 0;
  intTime = 0;
  DACsetting =0;
  for (int i = 0; i < 16; i ++){
    scalers[i] = 0;
  }
  return;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int vmeauxdata::Unpack(bankstructure bank){
  unsigned int* data = (unsigned int*)bank.Data;

  // Timing Board Data
  // data[0]=1  for version number (ignore for now)
  numAuxWords=(int)data[1];
  rampDelay = (int)data[2];
  intTime = (int)data[3];
  DACsetting= (int)data[4]; //DAC setting used for LED pulser scan
  // TIR data
  int mytirdata = data[5];
  dithering=(mytirdata & 0x100)>>8;     // channel 10? -- DOUBLE CHECK THIS
  helBit = (mytirdata & 0x10)>>4;
  cavpolBit = (mytirdata & 0x20)>>5;
  cavpowBit = (mytirdata & 0x40)>>6;
  mpsSignal = (mytirdata & 0x80)>>7;

  int point = numAuxWords+3;
  int point2 = point+data[point-1]+2;

  // Scaler data
  //if ((data[3]&0xcae56000) == 0xcae56000) {
  if(data[point] == 0xcae56000){
    for (int i = 0; i < 16; i ++){
      scalers[i] = data[point+i+1];
    }
  }
  else {              // Something has changed in the subbank
    cout << "\nERROR in vmeauxdata::Unpack: Scaler header does not appear where expected.";
  }

  if (data[point2] == 0xdae56000) {
    for (int i = 0; i < 16; i ++){
      IPscalers[i] = data[point2+1+i];
    }
  }
  else {              // Something has changed in the subbank
    cout << "\nERROR in vmeauxdata::Unpack: IP Scaler header does not appear where expected.";
  }


  //	DebugDump(2);

  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Which module information to dump?    0 = all
//                                      Multiples of 2 = Trigger Interface
//                                      Multiples of 3 = Timing Board
//                                      Multiples of 5 = Scalers
// e.g. 6 gives you both trigger interface and timing board
int vmeauxdata::DebugDump(int whichmodule){
  if (whichmodule%2 == 0){
    cout << "\nTIR data: helicity = " << helBit;
    cout << "\n          cavity power = " << cavpowBit;
    cout << "\n          cavity polarization direction = " << cavpolBit;
    cout << "\n          event length = " << evlen << endl;
    cout << "\n          mpsSignal bit = " << mpsSignal << endl;
    //		if (cavpow>1) cout << "\nHallelujah! Cavity was on!";
  }

  if (whichmodule%3 == 0){
    cout << "\nTiming Board data: ramp delay = " << rampDelay;
    cout << "\n                   integration time = " << intTime << endl;
  }

  if (whichmodule%5 == 0){

    for (int i = 0; i < 16; i++){
      cout << "\nScaler " << i << " count: " << scalers[i];
    }
    for (int i = 0; i < 16; i++){
      cout << "\nIPScaler " << i << " count: " << IPscalers[i];
    }
    cout << endl;
  }
  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
