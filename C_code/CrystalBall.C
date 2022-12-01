

TH1 *thing() {
  //TCutG *ne = (TCutG *)_file1->Get("ne");
  TH2D *h   = (TH2D *)_file0->Get("pid_addback_10ms_150");
  TH2D *hbg = (TH2D *)_file0->Get("pid_addback_10ms_150BG");

  TH1D *hx   = h->ProjectionX("hx",1,260,"[ne]");
  TH1D *hbgx = hbg->ProjectionX("hbgx",1,260,"[ne]");
  TH1D *chx = (TH1D *)hx->Clone("chx");
  chx->Add(hbgx,-1);
  chx->Rebin(4);
  
  return chx;
}

double Gaus(double *dim, double *par){
  // 3 parameters: a, x0, sigma;
  double a = par[0];
  double x0 = par[1];
  double sigma = par[2];
  double x = dim[0];  

  return a*TMath::Exp(-((x-x0)*(x-x0))/(2*sigma*sigma));

}

double DoubleGaus(double *dim, double *par){
  // 3 parameters: a, x0, sigma;
  double a1 = par[0];
  double x1 = par[1];
  double a2 = par[2];
  double x2 = par[3];
  double sigma = par[4];
  double x = dim[0];  

  double arr1[3] = {a1,x1,sigma};
  double arr2[3] = {a2,x2,sigma};
  return Gaus(dim,arr1) + Gaus(dim,arr2);

}

double CrystalBall(double *dim, double *par){
  double     A = par[0]; //par[0] = amplitude;
  double  mean = par[1]; //par[1] = mean; 
  double alpha = par[2]; //par[2] = alpha;
  double     n = par[3]; //par[3] = n;
  double sigma = par[4]; //par[4] = sigma;
  double     x = dim[0];

  if(sigma < 0.) return 0.;

  double z = (x-mean)/sigma;
  
  if(alpha < 0) z = -z;

  double abs_alpha = std::abs(alpha);
  
  if(z > -abs_alpha) return A*std::exp(-0.5*z*z); // Gaussian part;

  double nDivAlpha = n/abs_alpha;
  double AA = std::exp(-0.5*abs_alpha*abs_alpha);
  double B = nDivAlpha - abs_alpha;
  double arg = nDivAlpha/(B-z);  

  return A* AA * std::pow(arg,n);
}

double CrystalBall_BGAdd(double *dim, double *par){
  double     A = par[0]; //par[0] = amplitude;
  double  mean = par[1]; //par[1] = mean; 
  double alpha = par[2]; //par[2] = alpha;
  double     n = par[3]; //par[3] = n;
  double sigma = par[4]; //par[4] = sigma;
  double    bg = par[5]; //par[5] = background;
  double     x = dim[0];

  if(sigma < 0.) return 0.;

  double z = (x-mean)/sigma;
  
  if(alpha < 0) z = -z;

  double abs_alpha = std::abs(alpha);
  
  if(z > -abs_alpha) return A*std::exp(-0.5*z*z) + bg; // Gaussian part;

  double nDivAlpha = n/abs_alpha;
  double AA = std::exp(-0.5*abs_alpha*abs_alpha);
  double B = nDivAlpha - abs_alpha;
  double arg = nDivAlpha/(B-z);  

  return A* AA * std::pow(arg,n) + bg;
}


double DoubleCrystalBall(double *dim, double *par){
  double     A1 = par[0]; //par[0] = amplitude1;
  double  mean1 = par[1]; //par[1] = mean1; 
  double     A2 = par[2]; //par[2] = amplitude2;
  double  mean2 = par[3]; //par[3] = mean2; 
  double alpha  = par[4]; //par[4] = alpha;
  double     n  = par[5]; //par[5] = n;
  double sigma  = par[6]; //par[6] = sigma;
  double     x  = dim[0];

  double arr1[5] = {A1, mean1, alpha, n, sigma};
  double arr2[5] = {A2, mean2, alpha, n, sigma};

  return CrystalBall(dim,arr1) + CrystalBall(dim,arr2);

}


void quickcommand(TH1 *hist=0, double lower=10000, double upper=13000){
  //thing();
  if(hist == 0) hist = thing(); 
  TF1 *g  = new TF1("g",DoubleGaus,10000,13500,5);
  TF1 *g1 = new TF1("g1",Gaus,10000,13500, 3);
  TF1 *g2 = new TF1("g2",Gaus,10000,13500, 3);

  g->SetParName(0,"A1");
  g->SetParName(1,"C1");
  g->SetParName(2,"A2");
  g->SetParName(3,"C2");
  g->SetParName(4,"Sigma");

  //double upper = 12750;
  //double lower = 10000;
 
  g->SetParameters(200,12220,20,11970,350);
  g->SetParLimits(0,0,100000);
  g->SetParLimits(1,lower,upper);
  g->SetParLimits(2,0,100000);
  g->SetParLimits(3,lower,upper);
  g->SetParLimits(4,0,100000);

  g1->SetLineColor(kBlue);
  g2->SetLineColor(kGreen);

  hist->Fit(g,"","",lower,upper);
  g1->SetParameter(0,g->GetParameter(0));
  g1->SetParameter(1,g->GetParameter(1));
  g1->SetParameter(2,g->GetParameter(4));
  g2->SetParameter(0,g->GetParameter(2));
  g2->SetParameter(1,g->GetParameter(3));
  g2->SetParameter(2,g->GetParameter(4));
  

  hist->GetListOfFunctions()->Add(g1);
  hist->GetListOfFunctions()->Add(g2);
  gPad->Modified();
  gPad->Update();

}


void DoubleCrystalBallFit(TH1 *hist=0, double lower=10000, double upper=13500){
  if(hist == 0) hist = thing();
  TF1 *f  = new TF1("f", DoubleCrystalBall,lower,upper,7);
  TF1 *f1 = new TF1("f1",CrystalBall,lower,upper,5); // 30Ne
  TF1 *f2 = new TF1("f2",CrystalBall,lower,upper,5); // 31Ne
  
  f->SetParName(0,"A1");
  f->SetParName(1,"C1");
  f->SetParName(2,"A2");
  f->SetParName(3,"C2");
  f->SetParName(4,"alpha");
  f->SetParName(5,"n");
  f->SetParName(6,"sigma");

  f->SetParameters(300, 12200, 30, 11970, 1.6, 0.8, 180);
  f->SetParLimits(0,0,10000);
  f->SetParLimits(1,lower,upper);
  f->SetParLimits(2,0,1000);
  f->SetParLimits(3,lower,upper);
  f->SetParLimits(4,0,10000);
  f->SetParLimits(5,0,10000);
  f->SetParLimits(6,0,10000);
 
  f->FixParameter(4,1.6);
  f->FixParameter(5,0.78);
  f->FixParameter(6,179.2);
   
  hist->Fit(f,"","",lower,upper);
  f1->SetParameter(0,f->GetParameter(0));
  f1->SetParameter(1,f->GetParameter(1));
  f1->SetParameter(2,f->GetParameter(4));
  f1->SetParameter(3,f->GetParameter(5));
  f1->SetParameter(4,f->GetParameter(6));
  f2->SetParameter(0,f->GetParameter(2));
  f2->SetParameter(1,f->GetParameter(3));
  f2->SetParameter(2,f->GetParameter(4));
  f2->SetParameter(3,f->GetParameter(5));
  f2->SetParameter(4,f->GetParameter(6));
  f1->SetLineColor(kBlue);
  f2->SetLineColor(kGreen);
  
  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  gPad->Modified();
  gPad->Update();
}








