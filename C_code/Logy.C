


TGraph* Logy(TH1 *th1){
  
  TGraph *gr = 0;
  std::vector<double> myx;
  std::vector<double> myy;
  std::vector<double> mydx;
  std::vector<double> mydy;

  long entries = th1->GetNbinsX();
  long x = 1;

  for(x=1;x<=entries;x++){
    double y = th1->GetBinContent(x);
    if(y<=0) y = 1e-6;
    double logy = TMath::Log(y);
    if(logy<-1){
      printf("bin = %lu \t content = %f \t log = %f\n",x,logy,y);
      continue;
    }
    myx.push_back(th1->GetBinCenter(x));
    myy.push_back(logy);
    //mydx.push_back(0);                    
    //mydy.push_back(1/sqrt(y));           
  }
  gr = new TGraph(myx.size(),&myx[0],&myy[0]);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(kFullCircle);
  gr->Fit("pol1");

  return gr;

}


TH1D* logy(TH1 *th1){
  TH1D *cth1 = (TH1D *)th1->Clone("cth1");
  cth1->Reset();
  long entries = th1->GetNbinsX();
  for(long x=1;x<=entries;x++){
    double y = th1->GetBinContent(x);
    if(y<=0) y = 1e-6;
    y = TMath::Log(y);
    if(y<-3){
      printf("bin = %lu \t content = %f \t log = %f\n",x,TMath::Exp(y),y);
      continue;
    }
    cth1->SetBinContent(x,y);
  }
  //TF1 *fit = new TF1("fit","[0]-log(2.718)*(-0.693/[1])*x");
  //double max = cth1->GetBinContent(cth1->GetMaximumBin());
  //fit->SetParameters(max,5);
  //cth1->Fit(fit);
  cth1->Fit("pol1");
  return cth1;

}
