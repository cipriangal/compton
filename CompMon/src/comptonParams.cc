//  comptonParams.cc
//  Special header file used to collect parameters that are edited in this file
//
#include "comptonParams.h"
comptonParams::comptonParams(){
  return;
}
void comptonParams::newRun(int run) {
  // update parameters, including run-dependent parameters
  //
  runnum=run;
  channel_calorimeter = 0; // compton calorimeter FADC channel
  last_event = 1000000;
  sum_size = 3;
  ped_value = 2391.;  //changed for test data gbf

  pedestal_horizontalfinger = 3687;
  pedestal_verticalfinger = 3701;
  sthr1 = ped_value - 10.;
  sthr2 = ped_value - 60.;
  sum_pedchannel = ped_value * sum_size;
  sum_sthr1 = sthr1 * sum_size;
  pulse_sep = 120;
  pre_thr = 5;
  post_thr = 20;
  start_deadtime = 10;
  stop_deadtime = 10;
  snap_num = 200;
  sample_time = 0.000000005; 	// 5 ns per sample
  cavity_power_on_thresh = 900;		// Watts, Oct 13, 2010, for DVCS
  cavity_power_off_thresh = 500;		// Watts, Oct 13, 2010, for DVCS
  cavity_pol_on_thresh = 80; //Jan 18, 2010
  cavity_pol_off_thresh = 40; //Sep 23, 2009
  bcm_calibration_factor = 2.5*3*2.*7.68E-5*.912*1.29/10;	// microamps per normalized bcm scaler -- calibrated to epbcmu3, April 25, 2010
  cav_pow_calibration_factor = 0.0329;  //calibrate normalized cavity power to watts, May 7, 2010
  bcm_cutoff = 1;	// in microamps
  high_bcm_cutoff = 35;	// in microamps, set on April 2, 2010 for PREX running 

  HV_cutoff = 0; //set Mar 4, 2012, for g2p 

  helicity_delay=0;  //  set for debuggin
  helicity_structure=8;  //neede if we shift to a delayed helicity
  helicity_bits=30;
/*
  if(run < 22438 || (run>22618 && run<22641) || (run>22690 && run<22711) || (run>22796&& run<22802) || run>22804)helicity_delay = 8;         // for HAPPEX and PVDIS, moller, and 120Hz running
  else helicity_delay = 16;         // for PREX
  //helicity_delay = 8;        
  if(run < 20778) helicity_bits = 24;         // for HAPPEX
  else helicity_bits = 30;         // for PVDIS and PREX
  if(run < 22438 || (run>22618 && run<22641) || (run>22690 && run<22711) || (run>22796 && run<22802) || run>22804) helicity_structure = 4;   // for HAPPEX and PVDIS
  else helicity_structure = 8;   // for PREX
*/
return;
}

