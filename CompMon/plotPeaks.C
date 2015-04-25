void plotPeaks(){
  // read snapshot info from tree and plot'
  //makes histogram of pulse peaks (max deviation from pedestal)
  // also histos channels 1 and to for pedestal check
  //
  TFile *rootfile=new TFile("lnkCompMon.output");
  TTree *t1=(TTree*)rootfile->Get("snapshots");
  // Tree variables
  int numSamples;
  int mpsCoda;
  Float_t samples[1000];
  int laserState;
  int beamOn;
  int peakRAU;  //peak in Raw Adc Units
  t1->SetBranchAddress("numSamples",&numSamples);
  t1->SetBranchAddress("snap",&samples);
  t1->SetBranchAddress("mpsCoda",&mpsCoda);
  t1->SetBranchAddress("laserState",&laserState);
  t1->SetBranchAddress("beamOn",&beamOn);

  TH1F* hPeak= new TH1F("hPeak","Counts vs Wave Form  Peak (small RAU-> Big Peak",
			1024,0.,4096.);
  TH1F* hPed0= new TH1F("hPed0","WF chan 0",4096,0.,4096.);
  TH1F* hPed1= new TH1F("hPed1","WF chan 1",4096,0.,4096.);
  TH1F* hPed2= new TH1F("hPed2","WF chan 2",4096,0.,4096.);

  Int_t ntriggers=(Int_t)t1->GetEntries();
  printf("Number of waveform samples in file = %d\n",ntriggers);
  for(int i=0; i<ntriggers; i++){
    t1->GetEntry(i);

    peakRAU=1.E8;
    for(int j=0; j<numSamples; j++){
      if(samples[j]<peakRAU) peakRAU=samples[j];
    }
    hPeak->Fill(peakRAU);
    hPed0->Fill(samples[0]);
    hPed1->Fill(samples[1]);
    hPed2->Fill(samples[2]);
  }
  TCanvas* c1=new TCanvas("c1","WF Peaks");
  hPeak->Draw();
  int peakChannel=hPed1->GetMaximumBin();

  Double_t xmax=hPed1->GetBinCenter(peakChannel);
  printf("Chan 1 Max counts at x=%d\n",xmax);  
  hPed0->SetAxisRange(xmax-10.,xmax+10.);
  hPed1->SetAxisRange(xmax-10.,xmax+10.);
  hPed1->SetLineColor(2);
  hPed2->SetAxisRange(xmax-10.,xmax+10.);
  hPed2->SetLineColor(3);
  
  TCanvas* c2=new TCanvas("c2","FADC Pedestals", 20,20,800,300);
  hPed0->Draw();
  hPed1->Draw("same");
  hPed2->Draw("same");

}

