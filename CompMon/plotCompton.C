void plotCompton(){
  // read snapshot info from tree and plot'
  //makes histogram of pulse peaks (max deviation from pedestal)
  // also histos channels 1 and to for pedestal check
  //
  //TFile *rootfile=new TFile("lnkCompMon.output");
  TCanvas* c1=new TCanvas("c1","Triggered Pulse Data");


  hNorm_sums_asym->SetMinimum(-0.5);
  hNorm_sums_asym->SetMaximum(+0.5);
  c1->Divide(1,3);
  c1->cd(1);
  hTrig_sums_laserOn->Draw();
  c1->cd(2);
  hTrig_sums_laserOff->Draw();
  c1->cd(3);
  hNorm_sums_asym->Draw();
}

