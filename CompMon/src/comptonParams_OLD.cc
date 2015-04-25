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
  //cavity_power_on_thresh = 270;		// RT reading (VtoF)
  //cavity_power_on_thresh = 240;		// RT reading (VtoF)
  //cavity_power_on_thresh = 180;		// RT reading (VtoF) Sep 23, 2009
  if((runnum>19624&&runnum<19630) )
    {
      cavity_power_on_thresh = 60;		// Addition a couple of low power runs, Dec 16, 2009
    }
  /*
    else if(runnum > 22416 && runnum < 22438){
    cavity_power_on_thresh = 700;		// 30Hz run with new green laser, added Apr 13, 2010
    }
  */
  else if(run < 22438 ){
    cavity_power_on_thresh = 200;		// Watts, June 1, 2010, for HAPPEX and PVDIS
  }
  else if(run>23295) {
    cavity_power_on_thresh = 900;		// Watts, Oct 13, 2010, for DVCS
  }
  else{
    //cavity_power_on_thresh = 220;		// RT reading (VtoF) Nov 9, 2009
    cavity_power_on_thresh = 400;		// Watts, May 7, 2010
  }
  //cavity_power_off_thresh = 170;
  //cavity_power_off_thresh = 160;
  //cavity_power_off_thresh = 150; //Sep 23, 2009
  //cavity_power_off_thresh = 100; //Sep 23, 2009
  //cavity_power_off_thresh = 50; //Apr 28, 2010
  if(run < 22438)cavity_power_off_thresh = 160; //June 1, 2010, for HAPPEX and PREX
  else if(run>23295) {
    cavity_power_off_thresh = 500;		// Watts, Oct 13, 2010, for DVCS
  }
  else cavity_power_off_thresh = 250;//Watts, May 7, 2010
  //cavity_pol_on_thresh = 90; //Sep 23, 2009
  cavity_pol_on_thresh = 80; //Jan 18, 2010
  cavity_pol_off_thresh = 40; //Sep 23, 2009
  //bcm_calibration_factor = 2.*7.68E-5;	// microamps per normalized bcm scaler
  //bcm_calibration_factor = 3*2.*7.68E-5;	// microamps per normalized bcm scaler -- changed to U1 bcm and added factor of three, Sep 24, 2009
  if(run<19677){
    bcm_calibration_factor = 2.*7.68E-5*.912;	// microamps per normalized bcm scaler -- U3, Dec 14, 2009
  }
  else if (run>22436 && run<23352){
    bcm_calibration_factor = 3*2.*7.68E-5*.912*1.29;	// microamps per normalized bcm scaler -- calibrated to epbcmu3, April 25, 2010
  }
  else if (run>23351){
    bcm_calibration_factor = 2.5*3*2.*7.68E-5*.912*1.29/10;	// microamps per normalized bcm scaler -- calibrated to epbcmu3, April 25, 2010
  }
  else{
    bcm_calibration_factor = 3*2.*7.68E-5*.912;	// microamps per normalized bcm scaler -- calibrated to epbcmu3, Oct 11, 2009
  }
  cav_pow_calibration_factor = 0.0329;  //calibrate normalized cavity power to watts, May 7, 2010
  bcm_cutoff = 1;	// in microamps
  //high_bcm_cutoff = 85;	// in microamps, set on Nov 28, 2009 for HAPPEX analysis
  high_bcm_cutoff = 35;	// in microamps, set on April 2, 2010 for PREX running 
  if(run > 22961 && run < 22970) high_bcm_cutoff = 1.5;
  //high_bcm_cutoff = 4;	// in microamps, set on May 11, 2010 for low current PREX running 
  if(run>23295 && run < 24902) high_bcm_cutoff = 1; // set on Sep 29, 2010 for DVCS
  if(run > 24901) high_bcm_cutoff = -1; // set on Mar 3, 2012 for g2p
  //HV_cutoff = -1550; //set Feb 11, 2010, for HAPPEX
  //HV_cutoff = -1350; //set Apr 20, 2010, for early PREX
  HV_cutoff = 0; //set Mar 4, 2012, for g2p 
  //    helicity_delay = 0;       // for d2n
  if(run < 22438 || (run>22618 && run<22641) || (run>22690 && run<22711) || (run>22796&& run<22802) || run>22804)helicity_delay = 8;         // for HAPPEX and PVDIS, moller, and 120Hz running
  else helicity_delay = 16;         // for PREX
  //helicity_delay = 8;        
  if(run < 20778) helicity_bits = 24;         // for HAPPEX
  else helicity_bits = 30;         // for PVDIS and PREX
  if(run < 22438 || (run>22618 && run<22641) || (run>22690 && run<22711) || (run>22796 && run<22802) || run>22804) helicity_structure = 4;   // for HAPPEX and PVDIS
  else helicity_structure = 8;   // for PREX
return;
}

