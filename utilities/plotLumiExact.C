const double pi = acos(-1);
const double hplanck = 6.6261e-34; // J*s
const double clight = 299792458; //m/s
const double r0 = 2.81794e-15;//m
const double m2toBarn = 1e-28;
const double J2eV = 6.242e18;

// playtime
// double eSigmaX = 100e-6; //m;
// double eSigmaY = 100e-6; //m;should try sqrt(10-> emittance; 5-> beta*)
// double eSigmaL = 13e-12*clight; //m 
// double eFreq = 499e6; 
// double nElectron = 1/1.6e-19 * 65e-6;//A/eCharge
// double nBunches = 1; // 1320*78000 => bunch spacing of ~10ns


// double lambda = 532e-9; //m
// double gSigmaX = 100e-6; //m
// double gSigmaY = 100e-6; //m
// double gSigmaL = 12e-12 * clight; //3.6e-3 m -- 12ps FWHM; 4ps is 1.3mm
// double lPower = 10; //W
// double lBunchEnergy = lPower/(eFreq*nBunches);
// double nGammaBunch = lBunchEnergy * lambda / (hplanck * clight) ;//number;

//Calcs for eic
double eSigmaX = 400e-6; //m;
double eSigmaY = 400e-6/sqrt(50); //m;should try sqrt(10-> emittance; 5-> beta*)
double eSigmaL = 0.01; //m 
double eFreq = 78e3; 
double nElectron = 6.2e10; //number; 
double nBunches = 1320; // 1320*78000 => bunch spacing of ~10ns


double lambda = 532e-9; //m
double gSigmaX = 100e-6; //m
double gSigmaY = 100e-6; //m
double gSigmaL = 12e-12 * clight; //3.6e-3 m -- 12ps FWHM; 4ps is 1.3mm
double lPower = 10; //W
double lBunchEnergy = lPower/(eFreq*nBunches);
double nGammaBunch = lBunchEnergy * lambda / (hplanck * clight) ;//number;

double calcLcw(double*,double*);
double calcLpulsed(double*,double*);
double calcRatioPulseWave(double*,double*);

void plotLumiExact(){
  TCanvas *c1=new TCanvas("c1","c1");
  auto frame1 = c1->DrawFrame(0,1e24,10,1e38);
  frame1->GetYaxis()->SetTitle("Luminosity [1/cm^2 /s]");
  frame1->GetXaxis()->SetTitle("Crossing angle [deg]");

  cout<<"Pulse: energy per bunch: "<<lBunchEnergy<<endl;
  cout<<"Pulse: number of photons in a bunch: "<<nGammaBunch<<endl;

  TF1 *pulse = new TF1("pulse",calcLpulsed,0,10,0);
  pulse->SetLineColor(2);
  pulse->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);
  frame1->GetXaxis()->SetRangeUser(0,5);

}

//luminosity for a CW laser
double calcLcw(double *x, double *par){
  double xx = x[0]*pi / 180;
  double term1 = (1+cos(xx))/(std::sqrt(2*pi) * sin(xx) );
  double term2 = (nElectron * eFreq) * nGammaBunch*(eFreq*nBunches) / clight;
  double term3 = std::sqrt( (eSigmaX*eSigmaX + gSigmaX*gSigmaX) + (eSigmaY*eSigmaY + gSigmaY*gSigmaY) );

  return term1 * term2 / term3 /1e4; //1e4 to get to 1/cm^2/s; 
}

//luminosity for a pulsed laser
double calcLpulsed(double *x, double *par){
  double xx = x[0]*pi / 180 /2; // theta/2

  double term1 = eFreq * nElectron * nGammaBunch;
  double term4 = cos(xx)/(2*pi);
  double term2 = std::sqrt( eSigmaX*eSigmaX + gSigmaX*gSigmaX );
  double term3 = std::sqrt( cos(xx)*cos(xx) * (eSigmaY*eSigmaY + gSigmaY*gSigmaY) + sin(xx)*sin(xx) * (eSigmaL*eSigmaL + gSigmaL*gSigmaL) );

  return term1 * term4 / term2 / term3 /1e4; //1e4 to get to 1/cm^2/s; 
}

double calcRatioPulseWave(double *x, double *par){
  return calcLpulsed(x,par)/calcLcw(x,par);
}
