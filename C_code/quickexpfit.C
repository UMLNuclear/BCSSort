double MYExp(double *x, double *par){

  double lam1 = 0.693/par[1];
  return par[0]*TMath::Exp(-lam1*x[0]);

}

double Bateman(double *x, double *par, int n=1){
  double mul = 1;
  double sum = 0;
  for(int i=1;i<=n;i++){
    double lami = 0.693/par[i];
    if(i>=2 && (i<=n-1)) mul = mul*lami;
    double denominator = 1;
    for(int j=1;j<=n;j++){
      if(j==i) continue;
      double lamj = 0.693/par[j];
      denominator = denominator * (lamj - lami);
    }
    sum += TMath::Exp(-lami*x[0])/denominator; 
  }
  if(n==1) return par[0]*mul*sum;
  return (0.693/par[n])*(par[0]*mul*sum);
}


double MYFit(double *x, double *par) {
  double value = par[0] + MYExp(x,par+1) + MYExp(x,par+3) + MYExp(x,par+5);// + MYExp(x,par+7); + MYExp(x,par+9);

  return value;
}

double Draw_fit(double *x, double *par){
  double sum = 0;
  //if(par[0]==0) sum = par[1] + MYExp(x,par+2);
  sum = par[1] + MYExp(x,par+2);
  //else sum = par[1] + Bateman(x,par+2,par[0]);
  return sum;
}



TGraph *convert_th1_tgr(TH1D *th1, int binnum){

  TGraphErrors *gr = 0;
  std::vector<double> myx;
  std::vector<double> myy;
  std::vector<double> mydx;
  std::vector<double> mydy;
  for(int i=2;i<=binnum;i++){
    myx.push_back(th1->GetBinCenter(i));
    myy.push_back(th1->GetBinContent(i));
    mydx.push_back(0);                    // X-error
    mydy.push_back(sqrt(myy.back()));     // Y-error
  }
  gr = new TGraphErrors(myx.size(),
      &myx[0],
      &myy[0],
      &mydx[0],
      &mydy[0]);

  return gr;
}


TGraph* funfit(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2.){
  int npar = 7;
  
  std::vector<double> presult;
  TH1D *th1 = (TH1D *)hist1->Clone("th1");
  th1->Rebin(rebin);
  int nbins = th1->GetNbinsX();
  if(binnum<0) binnum = nbins;
  if(binnum>nbins) { printf("WARNING: Fitting range is beyond the histogram!!!\n"); }
  TGraph *gr = convert_th1_tgr(th1, binnum); 
  double Ymax = th1->GetMaximum();
  double Ymin = th1->GetMinimum();
  
  TF1 *fx = new TF1("fx",MYFit,0,binnum*rebin/10,npar);

  // Initila Pars Guess
  double p[npar];
  p[0] = 50000.;  
  p[1] = 10000;   
  p[2] = 300; 
  p[3] = 5000;   
  p[4] = 7.3;
  p[5] = 4500; 
  p[6] = 48;
  //p[7] = 500;
  //p[8] = 44.;
  //p[9] = 250;
  //p[10] = 337.;
  for(int i=0;i<npar;i++){
    fx->SetParameter(i,p[i]);
  }
  fx->SetParLimits(0,0,Ymax);
  fx->SetParLimits(1,0,Ymax);
  fx->SetParLimits(2,0,100000);
  fx->SetParLimits(3,0,Ymax);
  fx->SetParLimits(4,6.7,7.9);
  fx->SetParLimits(5,0,Ymax);
  fx->SetParLimits(6,46,50);
  //fx->SetParLimits(7,0,Ymax);
  //fx->SetParLimits(8,43,45);
  //fx->SetParLimits(9,0,Ymax);
  //fx->SetParLimits(10,310,360);

  gr->Fit(fx);
  presult.clear();
  for(int i=0;i<npar;i++){
    presult.push_back(fx->GetParameter(i));
  }
  double chi2 = fx->GetChisquare();
  double ndf = fx->GetNDF();  
  chi2 = chi2/ndf;

  while(chi2>1){
    for(int ipar=0;ipar<npar;ipar++){
      fx->SetParameter(ipar,presult[ipar]);
    }
    fx->SetParLimits(0,0,Ymax);
    fx->SetParLimits(1,0,Ymax);
    fx->SetParLimits(2,0,100000);
    fx->SetParLimits(3,0,Ymax);
    fx->SetParLimits(4,6.7,7.9);
    fx->SetParLimits(5,0,Ymax);
    fx->SetParLimits(6,46,50);
    //fx->SetParLimits(7,0,Ymax);
    //fx->SetParLimits(8,43,45);
    //fx->SetParLimits(9,0,Ymax);
    //fx->SetParLimits(10,310,360);

    gr->Fit(fx);


    double test_chi2 = fx->GetChisquare();
    test_chi2 = test_chi2/ndf;
    if((test_chi2 > chi2 && test_chi2<chi2_min) || (test_chi2<1)){
      break;  
    }
    if(test_chi2<chi2){
      presult.clear();
      for(int i=0;i<npar;i++){
        presult.push_back(fx->GetParameter(i));
      }
    }
    chi2 = test_chi2;
  }
  TF1 *f0 = new TF1("f0","[0]",0,binnum*rebin/10);     // constant bg
  f0->SetParameter(0,fx->GetParameter(0)+5e4);

  TF1 *f1 = new TF1("f1",Draw_fit,0,binnum*rebin/10,4);  // exp bg
  f1->SetParameter(0,1);                                 // return parent curve
  f1->SetParameter(1,presult[0]);                        // constant bg
  f1->SetParameter(2,presult[1]);                        // actvity
  f1->SetParameter(3,presult[2]);                        // t1/2

  TF1 *f2 = new TF1("f2",Draw_fit,0,binnum*rebin/10,4);  // 30Ne
  f2->SetParameter(0,1);                      // return parent curve
  f2->SetParameter(1,presult[0]+5e4);             // constant bg
  f2->SetParameter(2,presult[3]);             // actvity
  f2->SetParameter(3,presult[4]);             // parent t1/2

  //TF1 *f3 = new TF1("f3",Draw_fit,0,binnum*rebin/10,5);  // 30Na
  //f3->SetParameter(0,2);                      // model choose 
  //f3->SetParameter(1,p.first[0]);             // constant bg
  //f3->SetParameter(2,p.first[3]);             // acitivity
  //f3->SetParameter(3,p.first[4]);             // parent(30Ne) t1/2
  //f3->SetParameter(4,p.first[5]);             // daughter(30Na) t1/2

 
  f0->SetLineWidth(2);
  f0->SetTitle("flat bg");
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  f1->SetTitle("exp bg");
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(2);
  f2->SetTitle("30Ne");

  gr->GetListOfFunctions()->Add(f0);
  gr->GetListOfFunctions()->Add(f1);
  gr->GetListOfFunctions()->Add(f2);  
 
  //gr->GetYaxis()->SetRangeUser(presult[0]*0.995, Ymax*1.05);  

  return gr;

}
