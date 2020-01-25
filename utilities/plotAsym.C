const double pi = acos(-1);
const double hplanck = 4.135667696e-15; // eV*s
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

double detZposition = 60000; //mm
double transDetSize = 20; //mm

std::vector<TCanvas*> can;
TF1 *unpolXsec,*unpolXsecEwght,*al,*rg,*aWxSec,*aWxSecE,*a2WxSec;
TF1 *atUDrho,*atUDpos,*atRhoFixPhi;
TF2 *at,*atWxSec, *at2WxSec,*atWxSecE, *atAng, *atRhoY;
  
void updateConsts();
double calcUnpolXsection(double*,double*);
double calcUnpolXsectionE2wght(double*,double*);

double calcAL(double*,double*);
double calcALxSecWght(double*,double*);
double calcALxSecEWght(double*,double*);
double calcAL2xSecWght(double*,double*);

double calcAT(double*,double*);
double calcATrhoFixPhi(double*,double*);
double calcATrhoY(double*,double*);
double calcATang(double*,double*);
double calcAUDrho(double*,double*);
double calcAUDpos(double*,double*);
double calcATxSecWght(double*,double*);
double calcATxSecEWght(double*,double*);
double calcAT2xSecWght(double*,double*);

double rhoAng(double*,double*);
void drawStuff(int);

void plotAsym(){
  unpolXsec = new TF1("unpolXsec",calcUnpolXsection,0,1,1);
  unpolXsecEwght = new TF1("unpolXsecEwght",calcUnpolXsectionE2wght,0,1,1);
  can.push_back(new TCanvas("c1","c1"));
  can[0]->Divide(2);
  can[0]->cd(1);
  auto frame1 = gPad->DrawFrame(0,0,1,1);
  frame1->GetYaxis()->SetTitle("d#sigma/d#rho [barn]");
  frame1->GetXaxis()->SetTitle("#rho");
  can[0]->cd(2);
  auto frame1a = gPad->DrawFrame(0,0,1,1);
  frame1a->GetYaxis()->SetTitle("d#sigma/d#rho * #rho [barn]");
  frame1a->GetXaxis()->SetTitle("#rho");
  
  al=new TF1("al",calcAL,0,1,1);
  can.push_back(new TCanvas("c2","longitudinal Azz"));
  auto frame2 = can[1]->DrawFrame(0,-0.3,1,0.7);
  frame2->GetYaxis()->SetTitle("A_{long}");
  frame2->GetXaxis()->SetTitle("#rho");

  at=new TF2("at",calcAT,0,1,-pi,pi,1);
  can.push_back(new TCanvas("c3","transverse Azz"));
  can[2]->Divide(3);

  can.push_back(new TCanvas("c4","#rho dependence on angle"));
  auto frame3= can[3]->DrawFrame(0,1e-10,3.5,1e-2);
  rg = new TF1("rg",rhoAng,0,pi,2);

  can.push_back(new TCanvas("c5","weighted asymetries"));
  aWxSec = new TF1("aWxSec",calcALxSecWght,0,1,1);
  a2WxSec = new TF1("a2WxSec",calcAL2xSecWght,0,1,1);
  aWxSecE = new TF1("aWxSecE",calcALxSecEWght,0,1,1);
  can[4]->Divide(1,3);
  can[4]->cd(1);
  auto fr3 = gPad->DrawFrame(0,-0.1,1,0.35);
  fr3->SetTitle("x-section weighted asymmetry");
  fr3->GetXaxis()->SetTitle("#rho");
  can[4]->cd(2);
  auto fr4 = gPad->DrawFrame(0,0,1,0.15);
  fr4->SetTitle("x-section weighted asymmetry^2");
  fr4->GetXaxis()->SetTitle("#rho");
  can[4]->cd(3);
  auto fr5 = gPad->DrawFrame(0,-0.1,1,0.35);
  fr5->SetTitle("x-section*#rho weighted asymmetry");
  fr5->GetXaxis()->SetTitle("#rho");

  atAng=new TF2("atAng",calcATang,0,0.01,-pi,pi,2);
  can.push_back(new TCanvas("c6","transverse Azz"));
  can[5]->Divide(3);

  atUDrho=new TF1("atUDrho",calcAUDrho,0,1,1);
  can.push_back(new TCanvas("c7","transverse Azz UD rho"));
  auto fr6= can[6]->DrawFrame(0,0,1,0.3); 
  fr6->SetTitle("UD asymmetry");
  fr6->GetXaxis()->SetTitle("#rho");

  atUDpos=new TF1("atUDpos",calcAUDpos,-transDetSize,transDetSize,2);
  can.push_back(new TCanvas("c8","transverse Azz UD position"));
  auto fr7= can[7]->DrawFrame(-transDetSize,-0.04,transDetSize,0.04); 
  fr7->SetTitle(Form("UD asymmetry at z=%d m",int(detZposition/1000)));
  fr7->GetXaxis()->SetTitle("y position [mm]");

  atRhoFixPhi=new TF1("atRhoFixPhi",calcATrhoFixPhi,0,1,2);
  can.push_back(new TCanvas("c9","transverse Azz rho at #phi=0"));
  auto fr8= can[8]->DrawFrame(0,0,1,0.3); 
  fr8->SetTitle("AT asymmetry at #phi=0");
  fr8->GetXaxis()->SetTitle("#rho");

  atRhoY=new TF2("atRhoY",calcATrhoY,0,1,-transDetSize,transDetSize,2);
  can.push_back(new TCanvas("c10","transverse Azz rho, y"));
  can[9]->Divide(3);
  
  Ebeam = 5e9;
  updateConsts();
  drawStuff(0);

  Ebeam = 12e9;
  updateConsts();
  drawStuff(1);

  Ebeam = 18e9;
  updateConsts();
  drawStuff(2);

}

void drawStuff(int i){
  int color[3]={2,1,4};
  string tit[3]={"5 GeV","12 GeV","18 GeV"};

  can[0]->cd(1);
  unpolXsec->SetParameter(0,a);
  unpolXsec->SetLineColor(color[i]);
  unpolXsec->SetLineWidth(2);
  unpolXsec->DrawCopy("same");
  cout<<"Unpolarized xSection: "<<unpolXsec->Integral(0,1)<<" barn"<<endl;
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  can[0]->cd(2);
  unpolXsecEwght->SetParameter(0,a);
  unpolXsecEwght->SetLineColor(color[i]);
  unpolXsecEwght->SetLineWidth(2);
  unpolXsecEwght->DrawCopy("same");
  cout<<"Unpolarized xSection Ewght: "<<unpolXsecEwght->Integral(0,1)<<" barn"<<endl;
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[1]->cd();
  al->SetParameter(0,a);
  al->SetLineColor(color[i]);
  al->SetLineWidth(2);
  al->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[4]->cd(1);
  aWxSec->SetParameter(0,a);
  aWxSec->SetLineColor(color[i]);
  aWxSec->SetLineWidth(2);
  aWxSec->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  cout<<"Xsec weighted asym: "<<aWxSec->Integral(0,1)
      <<"\t <A>: "<<aWxSec->Integral(0,1)/unpolXsec->Integral(0,1)
      <<"\t <A>^2: "<<std::pow(aWxSec->Integral(0,1)/unpolXsec->Integral(0,1),2)<<endl;

  can[4]->cd(2);
  a2WxSec->SetParameter(0,a);
  a2WxSec->SetLineColor(color[i]);
  a2WxSec->SetLineWidth(2);
  a2WxSec->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  cout<<"Xsec weighted asym2: "<<a2WxSec->Integral(0,1)
      <<"\t <A^2>: "<<a2WxSec->Integral(0,1)/unpolXsec->Integral(0,1)<<endl;

  can[4]->cd(3);
  aWxSecE->SetParameter(0,a);
  aWxSecE->SetLineColor(color[i]);
  aWxSecE->SetLineWidth(2);
  aWxSecE->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  cout<<"Xsec*E weighted asym: "<<aWxSecE->Integral(0,1)
      <<"\t <EA>: "<<aWxSecE->Integral(0,1)/unpolXsecEwght->Integral(0,1)
      <<"\t <EA>^2: "<<std::pow(aWxSecE->Integral(0,1),2)/unpolXsecEwght->Integral(0,1)<<endl;

  
  can[2]->cd(i+1);
  at->SetParameter(0,a);
  at->SetTitle(Form("%s;#rho;#phi [rad]",tit[i].c_str()));
  at->DrawCopy("COLZ");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[5]->cd(i+1);
  atAng->SetParameters(a,Gamma);
  atAng->SetTitle(Form("%s;#theta [deg];#phi [rad]",tit[i].c_str()));
  atAng->DrawCopy("COLZ");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[3]->cd();
  rg->SetParameter(0,a);
  rg->SetParameter(1,Gamma);
  rg->SetLineColor(color[i]);
  rg->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLogy(1);

  can[6]->cd();
  atUDrho->SetParameter(0,a);
  atUDrho->SetLineColor(color[i]);
  atUDrho->SetLineWidth(2);
  atUDrho->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[8]->cd();
  atRhoFixPhi->SetParameters(a,0);
  atRhoFixPhi->SetLineColor(color[i]);
  atRhoFixPhi->SetLineWidth(2);
  atRhoFixPhi->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[7]->cd();
  atUDpos->SetParameters(a,Gamma);
  atUDpos->SetLineColor(color[i]);
  atUDpos->SetLineWidth(2);
  atUDpos->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  can[9]->cd(i+1);
  atRhoY->SetParameters(a,Gamma);
  atRhoY->SetTitle(Form("%s at z=%d m;#rho;y[mm] ",tit[i].c_str(),int(detZposition/1000)));
  atRhoY->DrawCopy("COLZ");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  
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

//unpolarized cross section for (x0)rho = E_gamma/E^max_gamma
double calcUnpolXsection(double *x, double *par){

  double term1 = std::pow(x[0],2) * std::pow(1 - par[0], 2) / ( 1 - x[0] * (1 - par[0]) );
  double term2 = std::pow( (1 - x[0] * (1+par[0])) / (1 - x[0] * (1-par[0])) , 2);

  return 2*pi*r0*r0*par[0] * (term1 + 1 + term2)/m2toBarn;
}

double calcUnpolXsectionE2wght(double *x, double *par){

  double term1 = std::pow(x[0],2) * std::pow(1 - par[0], 2) / ( 1 - x[0] * (1 - par[0]) );
  double term2 = std::pow( (1 - x[0] * (1+par[0])) / (1 - x[0] * (1-par[0])) , 2);

  return 2*pi*r0*r0*par[0] * (term1 + 1 + term2)/m2toBarn * x[0]*x[0];
}

/*********************************************************************************************/

//longitudinal asymmetry (Azz) for rho = E_gamma/E^max_gamma
double calcAL(double *x, double *par){
  double unpolXS = calcUnpolXsection(x,par)*m2toBarn;
  double term1 = 1 - x[0]*(1+par[0]);
  double term2 = 1 - 1/std::pow(1-x[0]*(1-par[0]),2);
  
  return 2*pi*r0*r0*par[0]/unpolXS * term1 * term2;
}

double calcALxSecWght(double *x, double *par){
  double asym = calcAL(x,par);
  double xSec = calcUnpolXsection(x,par);

  return asym*xSec;
}

double calcALxSecEWght(double *x, double *par){
  double asym = calcAL(x,par);
  double xSec = calcUnpolXsection(x,par);

  return asym*xSec*x[0];
}

double calcAL2xSecWght(double *x, double *par){
  double asym = std::pow(calcAL(x,par),2);
  double xSec = calcUnpolXsection(x,par);

  return asym*xSec;
}

/*********************************************************************************************/

//rho dependence on scattered photon angle
//x0=theta_photon[rad]; par0=a, par1=Gamma
double rhoAng(double *x, double *par){
  return 1/(1+par[0]*x[0]*x[0]*par[1]*par[1]);
}
//inverse of rhoAng; gives angle from rho
//x0=rho; par0=a, par1=Gamma
double angRho(double *x, double *par){
  return std::sqrt((1-x[0])/(x[0]*par[0]*par[1]*par[1]));
}

/*********************************************************************************************/

//transverse asymmetry (Azz) for phi (angle of the outgoing photon wrt e- polarization direction)
//function of rho(x0) and phi(x1); par0=a;
double calcAT(double *x, double *par){
  double unpolXS = calcUnpolXsection(x,par)*m2toBarn;

  double term1 = x[0]*(1-par[0])* std::sqrt(4*x[0]*par[0]*(1-x[0])) / (1-x[0]*(1-par[0]));

  return 2*pi*r0*r0*par[0]/unpolXS * cos(x[1]) * term1*100;
}

//function of rho(x0); par0=a, par1=phi;
double calcATrhoFixPhi(double *x, double *par){

  double xx[2] = {x[0],par[1]};
  double at = calcAT(xx,par);

  return at/100;//calcAT gives in percent
}


//function of theta_photon(x0[deg]) and phi(x1[rad]); par0=a, par1=Gamma
double calcATang(double *x, double *par){
  double xx = x[0]/180*pi;//convert to rad
  double rho = rhoAng(&xx,par);
  double xp[2] = {rho,x[1]};
  return calcAT(xp,par);
}


//function of rho(0-1),y(x1[mm]); par0=a, par1=Gamma
double calcATrhoY(double *x, double *par){

  double ang = angRho(x,par);
  double radius = tan(ang)*detZposition;
  double phi = acos(x[1]/radius);
  double xx[2] = {ang*180/pi, phi};
  double at = calcATang(xx,par);

  return at;
}

//function of rho(x0); par0=a;
double calcAUDrho(double *x, double *par){
  double vu(0),vd(0);
  int nu(0),nd(0);
  for(int i=0;i<=100;i++){
    double phi = -pi + i*2*pi/100;
    double xx[2] = {x[0],phi};
    double at = calcAT(xx,par);
    //double at = std::abs(calcAT(xx,par));
    if(phi < -pi/2 || phi>pi/2){ vd += at; nd ++;}
    else {vu += at; nu++;}
  }
  //cout<<vu<<"\t"<<vd<<endl;
  //std::cin.ignore();
  vu /= nu;
  vd /= nd;

  return (vu-vd)/200;//2 for mean, 100 calcAT gives res in percent
}

//function of yPosition(x0 [mm]); par0=a, par1=Gamma;
double calcAUDpos(double *x, double *par){
  double val(0);

  for(int i=0;i<=100;i++){
    double xPos = -transDetSize + i*2*transDetSize/100;
    double phi = atan2(xPos,x[0]);
    double ang = atan2(std::sqrt(x[0]*x[0]+xPos*xPos),detZposition);
    // double rho = rhoAng(&ang,par);
    // double unpolXS =  
    double xx[2] = {ang*180/pi,phi};
    double at = calcATang(xx,par);
    val += at;
  }
  val/=101;

  return val/100;//calcATang comes in percent
}

//yPos(x0[mm]), par0=a, par1=Gamma
double calcATxSecWght(double *x, double *par){
  double asym = calcAUDpos(x,par);
  // double rho = 
  double xSec = calcUnpolXsection(x,par);//rho, par0=a

  return asym*xSec;
}

double calcATxSecEWght(double *x, double *par){
  double asym = calcAT(x,par);
  double xSec = calcUnpolXsection(x,par);

  return asym*xSec*x[0];
}

double calcAT2xSecWght(double *x, double *par){
  double asym = std::pow(calcAT(x,par),2);
  double xSec = calcUnpolXsection(x,par);

  return asym*xSec;
}
