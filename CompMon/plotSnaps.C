void plotSnaps(){
  // read snapshot info from tree and plot'
  //setyo ti read riit fuke
  TFile *rootfile=new TFile("lnkCompMon.output");
  TTree *t1=(TTree*)rootfile->Get("snapshots");
  // Tree variables
  int numSamples;
  int mpsCoda;
  Float_t samples[1000];
  int laserState;
  int beamOn;

  t1->SetBranchAddress("numSamples",&numSamples);
  t1->SetBranchAddress("snap",&samples);
  t1->SetBranchAddress("mpsCoda",&mpsCoda);
  t1->SetBranchAddress("laserState",&laserState);
  t1->SetBranchAddress("beamOn",&beamOn);

  int gotoMPS=0;
  int skipMPS=0;
  int skipMPScounter=0;
  int oldMPS=0;

  //
  // setup for histogram
  TString cresponse;
  TCanvas* c1=new TCanvas("c1", "snapshots");
  TH1F* hsnap=new TH1F("hsnap","snapshot",500,0,500);
  Int_t ntriggers=(Int_t)t1->GetEntries();
  printf("Number of triggers in file = %d\n",ntriggers);
  for(int i=0; i<ntriggers; i++){
    t1->GetEntry(i);
    if(oldMPS!=mpsCoda){
      skipMPScounter++;
      oldMPS=mpsCoda;
    }
    if(skipMPScounter<skipMPS) continue;
    if(gotoMPS>0 && mpsCoda<gotoMPS) continue;
    skipMPScounter=0;
    gotoMPS=0;
    skipMPS=0;

    for(int j=0; j<numSamples; j++){
      hsnap->Fill(j,samples[j]);
    }


    //Get keyboard input
    hsnap->Draw();
    c1->Update();

    bool wait=true;
    printf("MPS  %6d Laser %2d BeamOn %2d >",mpsCoda, laserState, beamOn);
    while(wait){

      cin >> cresponse;
      cresponse.ToUpper();
      if(cresponse.Index("Q")==0) {
	return;
      }else if (cresponse.Index("N")==0) {
	wait=false;
      }else if(cresponse.Index("G")==0) {
	cin >>cresponse;   //get number MPS to Go To
	gotoMPS=cresponse.Atoi();
	wait=false;
      }else if(cresponse.Index("S")==0) {
	cin >>cresponse;   //get number MPS to Skip
	skipMPS=cresponse.Atoi();
	wait=false;
      }else if (cresponse.Index("?")==0){
	printf("q    Quit \n");
	printf("n    Next sample\n");
	printf("s x  Skip x MPSs\n");
	printf("g x  Goto MPS x\n");
	printf("?    Print commands\n");
	printf(" ?                           >");
      }else{
	printf(" ?                           >");
      }
    }

    hsnap->Reset();
  }
}
