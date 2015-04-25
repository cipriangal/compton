#include <stdio.h>
#include <iostream>
using namespace std;
#include <string>
#include "THaEpics.h"
#include "TString.h"

#ifndef epicscomptondata_h
#define epicscomptondata_h

class epicsComptonData
{
  private:
    float cavpow;
    float s1;
    float s2;
    float cavpoldir;
    float cavpolpercent;
    float hacbcm;
    float crystalHV;
    float crystalCurrent;
    float verticalFingerHV;
    float verticalFingerCurrent;
    float horizontalFingerHV;
    float horizontalFingerCurrent;
    float tablePosX;
    float tablePosY;
    float beamPosX;
    float beamPosY;
    float saclayRate;
    float beamEng;
    float transmit;
    int epics_counter;
    int ihwp_in;			// -1 = unknown, 0 = out, 1 = in
    int weinflip;
    float VWeinAngle;
    float Phi_FG;
    int helStruct;
    int helDelay;
    int helFreq;

  public:
    epicsComptonData(){
      cavpow = 0;
      cavpoldir = 0;
      cavpolpercent = 0;
      hacbcm = 0;
      s1 = 0;
      s2 = 0;
      crystalHV = -100.;		// hopefully this value will be obviously wrong
      crystalCurrent = -100.;	// as above
      verticalFingerHV = -100.;
      verticalFingerCurrent = -100.;
      horizontalFingerHV = -100.;
      horizontalFingerCurrent = -100.;
      tablePosX = -999.;
      tablePosY = -999.;
      beamPosX = -999.;
      beamPosY = -999.;
      saclayRate = -999.;
      beamEng = -999.;
      transmit = -999.;
      epics_counter = 0;
      ihwp_in = -1;			// unknown state
      weinflip = -1;
      VWeinAngle = 0;
      Phi_FG = 0;
      helStruct = -1;
      helDelay = -1;
      helFreq = -1;
    };

    float GetCavityPower()	{return cavpow;};
    float GetCavityPolarDirection()	{return cavpoldir;};
    float GetCavityPolarPercentage() {return cavpolpercent;};
    float GetBeamCurrent()	{return hacbcm;};
    float GetCrystalHV()	{return crystalHV;};	// in Volts
    float GetCrystalCurrent()	{return crystalCurrent;};		// in microAmps
    float GetVerticalFingerHV()  {return verticalFingerHV;};
    float GetVerticalFingerCurrent() {return verticalFingerCurrent;};
    float GetHorizontalFingerHV() {return horizontalFingerHV;};
    float GetHorizontalFingerCurrent()  {return horizontalFingerCurrent;};
    float GetTablePosX(){return tablePosX;};
    float GetTablePosY(){return tablePosY;};
    int GetIHWPState()		{return ihwp_in;};
    int GetWeinFlip()		{return weinflip;};
    int GetEPICSCounter()	{return epics_counter;};
    float GetS1Power()    {return s1;};
    float GetS2Power()    {return s2;};
    float GetBeamPositionX()  {return beamPosX;};
    float GetBeamPositionY()  {return beamPosY;};
    float GetSaclayRate() {return saclayRate;};
    float GetBeamEnergy() {return beamEng;};
    float GetTransmit()   {return transmit;};

    int Unpack(THaEpics *epics, uint32_t* codadata){
      TString spol;
      int evtype = codadata[1]>>16;
      int evnum = codadata[4];
      //  cout << "\nEvent #" << evnum << ": event type is " << evtype;
      if (evtype == 131)		// EPICS event
      {
        epics->LoadData(codadata,evnum);
        if (epics->IsLoaded("COMPTON_PW1PCAV_ca")){
          cavpow = epics->GetData("COMPTON_PW1PCAV_ca");
        }   // end cavpow loaded
        if (epics->IsLoaded("COMPTON_SU_POLAR_mo")){
          spol = epics->GetString("COMPTON_SU_POLAR_mo");
          cavpoldir = -1;
          if (spol.Strip(TString::kTrailing, ' ') == "RIGHT") cavpoldir = 1;
        }   // end cavpol loaded
        if (epics->IsLoaded("hac_bcm_average")){
          hacbcm = epics->GetData("hac_bcm_average");
        }   // end bcm loaded
        if (epics->IsLoaded("COMPTON_CAVPOLAR_ca")){
          cavpolpercent = epics->GetData("COMPTON_CAVPOLAR_ca");
        }
        if (epics->IsLoaded("COMPTON_PW1R_S1_ca")){
          s1 = epics->GetData("COMPTON_PW1R_S1_ca");
        }
        if (epics->IsLoaded("COMPTON_PW1R_S2_ca")){
          s2 = epics->GetData("COMPTON_PW1R_S2_ca");
        }
        if (epics->IsLoaded("MF2b1461_11ch1property.F")){
          crystalHV = epics->GetData("MF2b1461_11ch1property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch1property.E")){
          crystalCurrent = epics->GetData("MF2b1461_11ch1property.E");	
        }
        if (epics->IsLoaded("MF2b1461_11ch8property.F")){
          verticalFingerHV = epics->GetData("MF2b1461_11ch8property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch8property.E")){
          verticalFingerCurrent = epics->GetData("MF2b1461_11ch8property.E");
        }
        if (epics->IsLoaded("MF2b1461_11ch7property.F")){
          horizontalFingerHV = epics->GetData("MF2b1461_11ch7property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch7property.E")){
          horizontalFingerCurrent = epics->GetData("MF2b1461_11ch7property.E");
        }
        if (epics->IsLoaded("COM_DETPH_XCPOSai")){
          tablePosX = epics->GetData("COM_DETPH_XCPOSai");
        }
        if (epics->IsLoaded("COM_DETPH_YCPOSai")){
          tablePosY = epics->GetData("COM_DETPH_YCPOSai");
        }
        if (epics->IsLoaded("IPM1P02B.XPOS")){
          beamPosX = epics->GetData("IPM1P02B.XPOS");
        }
        if (epics->IsLoaded("IPM1P02B.YPOS")){
          beamPosY = epics->GetData("IPM1P02B.YPOS");
        } 
        if (epics->IsLoaded("compton:RATE_G")){
          saclayRate = epics->GetData("compton:RATE_G");
        }
        if (epics->IsLoaded("HALLA:p")){
          beamEng = epics->GetData("HALLA:p");
        }
        if (epics->IsLoaded("COMPTON_TRANSMIT_ao")){
          transmit = epics->GetData("COMPTON_TRANSMIT_ao");
        }
        if (epics->IsLoaded("HELPATTERNd")){
          helStruct = epics->GetData("HELPATTERNd");
        }
        if (epics->IsLoaded("HELDELAYd")){
          helDelay = epics->GetData("HELDELAYd");
        }
        if (epics->IsLoaded("HELFREQ")){
          helFreq = epics->GetData("HELFREQ");
        }
        if (epics->IsLoaded("IGL1I00OD16_16")){
          spol = epics->GetString("IGL1I00OD16_16");
          if (spol.Strip(TString::kTrailing, ' ') == "IN"){
            ihwp_in = 1;
          }
          else if (spol.Strip(TString::kTrailing, ' ') == "OUT"){ 
            ihwp_in = 0;
          }
          else{
            ihwp_in = -1;		// Unknown!
          }
        }
        if (epics->IsLoaded("VWienAngle") && epics->IsLoaded("Phi_FG")){
          VWeinAngle = epics->GetData("VWienAngle");
          Phi_FG = epics->GetData("Phi_FG");
          if((VWeinAngle < 100 && VWeinAngle > 80) && (Phi_FG < 100 && Phi_FG > 80)){
            weinflip = 1;
          }
          else if((VWeinAngle > -100 && VWeinAngle < -80 && Phi_FG < 100 && Phi_FG > 80) || (Phi_FG > -100 && Phi_FG < -80 && VWeinAngle < 100 && VWeinAngle > 80)){
            weinflip = 0;
          }
          else{
            weinflip = -1;
          }
        }
        epics_counter++;

        //		DebugDump();
      }      // end evtype == 131
      return 0;
    };
    
    int DebugDump(){
      cout << "\n----------EPICS DATA DUMP-----------";
      cout << "\nCavity power: " << cavpow;
      cout << "\nCavity polarization direction: " << cavpoldir;
      cout << "\nCavity polarization percentage: " << cavpolpercent;
      cout << "\nHigh voltage on crystal PMT: " << crystalHV;
      cout << "\nCurrent on crystal PMT: " << crystalCurrent;
      cout << "\nAverage beam current: " << hacbcm;
      cout << "\nInsertable half-waveplate state: " << ihwp_in;
      cout << endl;
       
      return 0;
    };
};

#endif			//epicsdomptondata_h
