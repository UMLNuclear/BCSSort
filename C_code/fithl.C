
TGraph *fithl(TH1D *th1, int rebin=8, int binnum=50){
  TGraphErrors *gr = 0;
  //TF1 fx("fx","[0]*TMath::Exp(TMath::Log(2)/[1]*x)+[2]");
  TF1 *fx = new TF1("fx","[0]*exp(-(0.693/[1])*x)+[2]",0,500);
  //TF1 *fx = new TF1("fx","[0]*exp(-(0.693/[1])*x)+[2]*exp(-(0.693)/[3]*x)+[4]*exp(-(0.693)/[5]*x)",0,300);
  //TF1 *fx = new TF1("fx","[0]*exp(-(0.693/[1])*x)+[2]*exp(-(0.693)/[3]*x)+[4]*exp(-(0.693)/[5])+[6]",0,300);
  TH1D *cth1 = (TH1D *)th1->Clone("cth1");
  cth1->Rebin(rebin);
  std::vector<double> myx;
  std::vector<double> myy;
  std::vector<double> mydx;
  std::vector<double> mydy;
  for(int i=1;i<=binnum;i++){
    myx.push_back(cth1->GetBinCenter(i));
    myy.push_back(cth1->GetBinContent(i));
    mydx.push_back(((double)rebin)*0.5*(1/3));
    mydy.push_back(sqrt(myy.back()));
    std::cout<<myx.back()<<"\t"
             <<myy.back()<<"\t"
             <<mydx.back()<<"\t"
             <<mydy.back()<<std::endl;
  } 
  gr = new TGraphErrors(myx.size(),
                        &myx[0],
                        &myy[0],
                        &mydx[0],
                        &mydy[0]);
  
  //fx->SetParLimits(6,0,500);
  fx->SetParameters(myy[0],8,myy.back());
  //fx->SetParameters(myy[0]/2,8,myy[0]/2,13,myy.back(),7000);
  //fx->SetParameters(myy[0]/10,8,myy[0],10,myy[0],15,myy.back());
  //fx->SetParameters(myy[0],8,myy[0],10,myy[0],10,myy[0],10,myy.back());
  
  //cout << fx->GetParameter(0) << "\t";
  //cout << fx->GetParameter(1) << "\t";
  //cout << fx->GetParameter(2) << endl;
  //fx->Print("all");
  //gr->Print();
  gr->Fit(fx);

  return gr;
}
