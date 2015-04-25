//  comptonParams.h
//  comptonParams deals with run-dependent parameters
//
#include <stdlib.h>
#include <TROOT.h>
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

#ifndef comptonParams_h
#define comptonParams_h

#define MAX_FADC_CHANNELS 8

class comptonParams
{
  private:
  int channel_calorimeter;          //fadc channel for calorimeter
  int channel_verticalfinger;      //fadc channel for vertical finger
  int channel_horizontalfinger;    //    for horizontal figure
    double ped_value;			// pedestal value parameter for calorimeter
    double sum_pedchannel;
    double pedestal_verticalfinger;
    double pedestal_horizontalfinger;
    double sthr1;				// near-pedestal sampling threshold
    double sthr2;				// far-pedestal sampling threshold
    int pulse_sep;				// enforced sample separation between pulses
    int pre_thr;				// number of samples to check before a pulse
    int post_thr;				// number of samples to check after a pulse
    int start_deadtime;			// deadtime at the start of a CODA event
    int stop_deadtime;			// deadtime at the end of a CODA event
    int channel;				// which FADC channel?
    int helicity_delay;     // delayed by how many MPSes?
    int helicity_bits;     // helicity (board) predictor bits
    int helicity_structure;     // quartet or octet?
    double sample_time;
    int last_event; 
    int sum_size;
    int snap_num;
    float bcm_calibration_factor;
    unsigned int bcm_cutoff;
    unsigned int high_bcm_cutoff;
    float HV_cutoff;

    float sum_sthr1;
    unsigned int cavity_power_on_thresh;	// How much power should the cavity have to be "on"?
    unsigned int cavity_power_off_thresh;
    float cavity_pol_on_thresh;	
    float cavity_pol_off_thresh;
    int runnum;
    float cav_pow_calibration_factor;


  public:
    comptonParams();  //constructor
    void newRun(int run); //update parameters  (possiblyi run-dependent)
    // Get parameters
    int GetCalorimeterChannel() {return channel_calorimeter;};
    double GetPedestal() {return ped_value;};
    double GetPedTop()  {return pedestal_verticalfinger;};
    double GetPedSide() {return pedestal_horizontalfinger;};
    double GetSamplingNearThreshold() {return sthr1;};
    double GetSamplingFarThreshold() {return sthr2;};
    int GetPulseSeparation() {return pulse_sep;};
    int GetPreThreshold() {return pre_thr;};
    int GetPostThreshold() {return post_thr;};
    int GetStartDeadtime() {return start_deadtime;};
    int GetStopDeadtime() {return stop_deadtime;};
    int GetFADCChannel() {return channel;};
    double GetSampleTime() {return sample_time;};
    int GetLastEvent() {return last_event;};
    int GetSumSize() {return sum_size;};
    int GetSnapNum() {return snap_num;};
    double GetSumPedChannel() {return sum_pedchannel;};
    float GetSumNearThreshold() {return sum_sthr1;};
    unsigned int GetCavPowOn()   {return cavity_power_on_thresh;};
    unsigned int GetCavPowOff()  {return cavity_power_off_thresh;};
    float GetCavPolOn()   {return cavity_pol_on_thresh;};
    float GetCavPolOff()  {return cavity_pol_off_thresh;};
    float GetBCMCalibration() {return bcm_calibration_factor;};
    //float GetU3BCMCalibration() {return U3_bcm_calibration_factor;};
    unsigned int GetBCMCutoff() {return bcm_cutoff;};
    unsigned int GetHighBCMCutoff() {return high_bcm_cutoff;};
    float GetHVCutoff() {return HV_cutoff;};
    int GetHelicityDelay()  {return helicity_delay;};
    int GetHelicityBits()  {return helicity_bits;};
    int GetHelicityStructure()  {return helicity_structure;};
    int GetRun()  {return runnum;};
    float GetCavPowCalibration() {return cav_pow_calibration_factor;};
    void SetTopPed(double x)  {pedestal_verticalfinger = x;};
    void SetSidePed(double x) {pedestal_horizontalfinger = x;};
    void SetCavPowOn(unsigned int x) {cavity_power_on_thresh = x;};
    void SetCavPowOff(unsigned int x)  {cavity_power_off_thresh = x;};
    void SetSamplingFarThreshold(double x) {sthr2 = ped_value - x;};
    void SetPulseSeparation(int x) {pulse_sep = x;};
    void SetPreThreshold(int x) {pre_thr = x;};
    void SetPostThreshold(int x) {post_thr = x;};
    void SetStartDeadtime(int x) {start_deadtime = x;};
    void SetStopDeadtime(int x) {stop_deadtime = x;};
    void SetFADCChannel(int x) {if (x < MAX_FADC_CHANNELS){channel = x;}};
    void SetSampleTime(double x) {sample_time = x;};
    void SetLastEvent(int x) {last_event = x;};
    //void SetSumSize(int x) {sum_size = x;
    // sum_pedchannel = ped_channel * sum_size;
    //sum_sthr1 = sthr1 * sum_size;};
    void SetSnapNum(int x) {snap_num = x;};
    void SetHelicityDelay(int x)  {helicity_delay = x;};
};
#endif
