const double pi = acos(-1);
//const double hplanck = 4.135667696e-15; // eV*s
const double hplanck = 6.6261e-34; // J*s
const double clight = 299792458; //m/s
const double r0 = 2.81794e-15;//m
const double m2toBarn = 1e-28;

//double lambda = 323e-9; //m
double lambda = 532e-9; //m
double Ebeam = 1e9; //eV
double melectron = 0.511e6; //eV/c^2

double Elaser= hplanck * clight / lambda; //eV
double Gamma = Ebeam / melectron;
double a = 1 / ( 1 + 4*Gamma*Elaser/melectron);//
double EgammaMax = 4*a*Elaser*std::pow(Gamma,2);//eV


// //this recoveres DGs plot (with a factor difference not sure which)
// // c1 range 1e28 to 1e31
// double eSigmaT = 100e-6; //m;
// double gSigmaT = 100e-6; //m
// double lPower = 1e3; //W
// double nElectron = 1/1.6e-19 * 50e-6; //A/e-charge; #/s
// double gSigmaL = 6e-12 * clight; //m -- 12ps FWHM
// double eSigmaL = 0.5e-12 *clight; //5ps*c
// double eFreq = 499e6; //499e6Hz CEBAF;


// // c1 range 1e24 to 1e27
double eSigmaT = 400e-6; //m; eRHIC 400e-6; CEBAF 100e-6
double gSigmaT = 100e-6; //m
double gSigmaL = 12e-12 * clight; //m -- 12ps FWHM
double eSigmaL = 13e-12 *clight; //m ~13 ps eRHIC; 0.5ps CEBAF;
double lPower = 1e3; //W
double eFreq = 78e3; //98e6Hz (*6 buckets) -eRHIC; 499e6Hz CEBAF;
double nElectron = 1/1.6e-19 * 10e-9; //A/e-charge; #/s (Elke used 7e9 e-, 10nC is 6e10)

double calcLcw(double*,double*);
double calcLpulsed(double*,double*);
double calcRatioPulseWave(double*,double*);

void plotLumi(){
  TCanvas *c1=new TCanvas("c1","c1");
  auto frame1 = c1->DrawFrame(0,1e24,10,1e38);
  frame1->GetYaxis()->SetTitle("Luminosity [1/cm^2 /s]");
  frame1->GetXaxis()->SetTitle("Crossing angle [deg]");

  TF1 *wave = new TF1("wave",calcLcw,0,10,0);
  wave->SetLineColor(4);
  wave->SetLineWidth(2);  

  TF1 *pulse = new TF1("pulse",calcLpulsed,0,10,0);
  pulse->SetLineColor(2);
  wave->SetLineWidth(2);  

  pulse->DrawCopy("same");
  wave->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);
  frame1->GetXaxis()->SetRangeUser(0.5,3);

  auto *c2=new TCanvas("c2","ratio");
  TF1 *rpw = new TF1("rpw",calcRatioPulseWave,0,10,0);
  rpw->DrawCopy();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  // TF1 *al=new TF1("al",calcAL,0,1,1);
  // TCanvas *c2=new TCanvas("c2","longitudinal Azz");
  // auto frame2 = c2->DrawFrame(0,-0.3,1,0.7);
  // frame2->GetYaxis()->SetTitle("A_{long}");
  // frame2->GetXaxis()->SetTitle("#rho");

}

//luminosity for a CW laser
double calcLcw(double *x, double *par){
  double xx = x[0]*pi / 180;
  double term1 = (1+cos(xx))/(std::sqrt(2*pi) * sin(xx) );
  double term2 = (nElectron * eFreq) * lPower * lambda / (hplanck * clight * clight);
  double term3 = std::sqrt(eSigmaT * eSigmaT + gSigmaT * gSigmaT);

  return term1 * term2 / term3 /1e4; //1e4 to get to 1/cm^2/s; 
}

//luminosity for a pulsed laser
double calcLpulsed(double *x, double *par){
  double xx = x[0]*pi / 180;
  double term1 = (1+cos(xx))/(2*pi * sin(xx) );
  double term2 = nElectron * lPower * lambda / (hplanck * clight);
  double term3 = std::sqrt(eSigmaT * eSigmaT + gSigmaT * gSigmaT);
  double term4 = std::sqrt(eSigmaL*eSigmaL + gSigmaL*gSigmaL + term3*term3 / std::pow(sin(xx/2),2));
  
  return term1 * term2 / term3 / term4 /1e4; //1e4 to get to 1/cm^2/s; 
}

double calcRatioPulseWave(double *x, double *par){
  return calcLpulsed(x,par)/calcLcw(x,par);
}
