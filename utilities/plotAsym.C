
const double pi = acos(-1);
const double hplanck = 4.135667696e-15; // eV*s
const double clight = 299792458; //m/s
const double r0 = 2.81794e-15;//m
const double m2toBarn = 1e-28;

double lambda = 532e-9; //m
double Ebeam = 1e9; //eV
double melectron = 0.511e6; //eV/c^2

double Elaser= hplanck * clight / lambda; //eV
double Gamma = Ebeam / melectron;
double a = 1 / ( 1 + 4*Gamma*Elaser/melectron);//
double EgammaMax = 4*a*Elaser*std::pow(Gamma,2);//eV

void updateConsts();
double calcUnpolXsection(double*,double*);
double calcAL(double*,double*);
double calcAT(double*,double*);

void plotAsym(){
  TF1 *unpolXsec = new TF1("unpolXsec",calcUnpolXsection,0,1,1);
  TCanvas *c1=new TCanvas("c1","c1");
  auto frame1 = c1->DrawFrame(0,0,1,1);
  frame1->GetYaxis()->SetTitle("d#sigma/d#rho [barn]");
  frame1->GetXaxis()->SetTitle("#rho");

  TF1 *al=new TF1("al",calcAL,0,1,1);
  TCanvas *c2=new TCanvas("c2","longitudinal Azz");
  auto frame2 = c2->DrawFrame(0,-0.3,1,0.7);
  frame2->GetYaxis()->SetTitle("A_{long}");
  frame2->GetXaxis()->SetTitle("#rho");

  updateConsts();
  c1->cd();
  unpolXsec->SetParameter(0,a);
  unpolXsec->SetLineWidth(2);
  unpolXsec->DrawCopy("same");
  c2->cd();
  al->SetParameter(0,a);
  al->SetLineWidth(2);
  al->DrawCopy("same");

  Ebeam = 11e9;
  updateConsts();
  c1->cd();
  unpolXsec->SetParameter(0,a);
  unpolXsec->SetLineColor(1);
  unpolXsec->SetLineWidth(2);
  unpolXsec->DrawCopy("same");
  c2->cd();
  al->SetParameter(0,a);
  al->SetLineColor(1);
  al->SetLineWidth(2);
  al->DrawCopy("same");
  
  Ebeam = 27e9;
  updateConsts();
  c1->cd();
  unpolXsec->SetParameter(0,a);
  unpolXsec->SetLineColor(4);
  unpolXsec->SetLineWidth(2);
  unpolXsec->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  c2->cd();
  al->SetParameter(0,a);
  al->SetLineColor(4);
  al->SetLineWidth(2);
  al->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  

  TF2 *at=new TF2("at",calcAT,0,1,-pi,pi,1);
  TCanvas *c3=new TCanvas("c3","transverse Azz");
  c3->Divide(3);

  c3->cd(1);
  Ebeam = 1e9;
  updateConsts();
  at->SetParameter(0,a);
  at->DrawCopy("COLZ");
  
  c3->cd(2);
  Ebeam = 11e9;
  updateConsts();
  at->SetParameter(0,a);
  at->DrawCopy("COLZ");
  
  c3->cd(3);
  Ebeam = 27e9;
  updateConsts();
  at->SetParameter(0,a);
  at->DrawCopy("COLZ");
  
}

void updateConsts(){
  Elaser= hplanck * clight / lambda; //eV
  Gamma = Ebeam / melectron;
  a = 1 / ( 1 + 4*Gamma*Elaser/melectron);
  EgammaMax = 4*a*Elaser*std::pow(Gamma,2);

  cout<<endl<<endl;
  cout<<"E laser: "<<Elaser<<" eV"<<endl;
  cout<<"E beam : "<<Ebeam <<" eV"<<endl;
  cout<<"gamma: "<<Gamma<<" "<<endl;
  cout<<"E gamma max: "<<EgammaMax<<" eV"<<endl;
  cout<<"a: "<<a<<" "<<endl;
}

//unpolarized cross section for rho = E_gamma/E^max_gamma
double calcUnpolXsection(double *x, double *par){

  double term1 = std::pow(x[0],2) * std::pow(1 - par[0], 2) / ( 1 - x[0] * (1 - par[0]) );
  double term2 = std::pow( (1 - x[0] * (1+par[0])) / (1 - x[0] * (1-par[0])) , 2);

  return 2*pi*r0*r0*par[0] * (term1 + 1 + term2)/m2toBarn;
}

//longitudinal asymmetry (Azz) for rho = E_gamma/E^max_gamma
double calcAL(double *x, double *par){
  double unpolXS = calcUnpolXsection(x,par)*m2toBarn;
  double term1 = 1 - x[0]*(1+par[0]);
  double term2 = 1 - 1/std::pow(1-x[0]*(1-par[0]),2);
  
  return 2*pi*r0*r0*par[0]/unpolXS * term1 * term2;
}

//transverse asymmetry (Azz) for phi (angle of the outgoing photon wrt e- polarization direction)
//function of rho and phi
double calcAT(double *x, double *par){
  double unpolXS = calcUnpolXsection(x,par)*m2toBarn;

  double term1 = x[0]*(1-par[0])* std::sqrt(4*x[0]*par[0]*(1-x[0])) / (1-x[0]*(1-par[0]));

  return 2*pi*r0*r0*par[0]/unpolXS * cos(x[1]) * term1;
}
