

#include<cstdio>
#include<TGraph.h>
#include<TGraphErrors.h>






//================================= Fitting Model ==================================//
// simple exponential function;
// par[0]: constant scaling factor
// par[1]: tau (halflife)
// dim[0]: x value - time 
double MYExp(double *x, double *par){

  double lam1 = 0.693/par[1];
  return par[0]*TMath::Exp(-lam1*x[0]);

}

// Calculate activity of one decay = lambda *N(t);
// N(t) is calculated by Bateman equation;
// n = isotope number(parent = 1, daughter = 2, granddaughter = 3);
// when n = 1, return exponential function par[0]*exp(-lam1*x[0])
// par[0] = particle[n] activity;par[i] = half-life of particle_i(i>0)
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


// sum of activities for a decay chain;
// n =  # of isotopes decay in this decay chain;
// # of par = n+1(n half-lifes + 1 parent activity);
// par[0] = the only one activity for each chain
double sum_Bateman(double *x, double *par, int n){
  double sum = 0;
  for(int i=1;i<=n;i++){
    sum += Bateman(x,par,i);
  }
  return sum;
}

//Eventual Fitting Model
//Example: 2 decay paths: 32Na->32Mg->32Al & 32Na->31Mg
//two bg: constant + exp bg
double MYFit(double *x, double *par) {
  //par[0]: A_30Ne
  //par[1]: t1/2_30Ne
  //par[2]: t1/2_30Na
  //par[3]: t1/2_29Na
  //par[4]: t1/2_30Mg
  //par[5]: br_(30Ne_bn)
  //par[6]: br_(30Na_bn)
  //par[7]: A_exp bg  
  //par[8]: t1/2_exp bg
  //par[9]: br(30Ne_b2n)  

  double par1[2]; //30Ne
  double par2[3]; //30Na
  double par3[3]; //29Na
  double par4[4]; //30Mg

  par1[0]=par[0];
  par1[1]=par[1];
  
  par2[0]=par[0];
  par2[1]=par[1];
  par2[2]=par[2];

  par3[0]=par[0];
  par3[1]=par[1];
  par3[2]=par[3];
 
  par4[0]=par[0];
  par4[1]=par[1];
  par4[2]=par[2];
  par4[3]=par[4];

  double bg = par[7]*TMath::Exp(-0.693/par[8]*x[0]);
  
  double value = bg + Bateman(x, par1,1) + (1-par[5])*Bateman(x,par2,2) + par[5]*Bateman(x,par3,2) + (1-par[6])*(1-par[5])*Bateman(x,par4,3);
  
  return value;
}
//================================== Math model for fitting plot =================================//
//par[0] = 0 => bg plot;
//par[0] >=0 => the value of par[0] = isotope number (parent=1, daughter=2, granddaughter=3)
double Draw_fit(double *x, double *par){
  double sum = 0;

  //if(par[0]==0) sum = par[1]*(1 + par[2]*x[0] + par[3]*x[0]*x[0]);
  if(par[0]==0) sum = par[1]*(0+TMath::Exp(-0.693/par[2]*x[0]));
  else if(par[0]<0) sum = par[1]*TMath::Exp(-0.693/par[2]*x[0]) + par[3]*Bateman(x,par+4,fabs(par[0]));
  else sum = par[1] + par[2]*Bateman(x,par+3,par[0]);
  
  return sum;
}



//================================== Convert TH1D to TGraphErrors =========================================================//

TGraph *convert_th1_tgr(TH1D *th1, int binnum){

  std::ofstream ofile1;  
  ofile1.open("decaycurve.dat");
  
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
    ofile1<<myx.back()<<"\t"
      <<myy.back()/100000.<<"\t"
      <<mydx.back()<<"\t"
      <<mydy.back()/100000.<<std::endl;
  }
  gr = new TGraphErrors(myx.size(),
      &myx[0],
      &myy[0],
      &mydx[0],
      &mydy[0]);

  return gr;
}

TGraph *fithl(TH1D *hist1, int rebin=1, int binnum=-1){


  TH1D *th1 = (TH1D *)hist1->Clone("th1");
  th1->Rebin(rebin);
  double Ymax = th1->GetMaximum();
  int nbins = th1->GetNbinsX();
  if(binnum<0) binnum = nbins;
  if(binnum>nbins) { printf("WARNING: Fitting range is beyond the histogram!!!\n"); }

  int npar = 9;
  TGraph *gr = convert_th1_tgr(th1, binnum); 
  TF1 *fx = new TF1("fx",MYFit,0,300,npar);

  fx->SetParameter(0,13915.417255);                                                 
  fx->SetParameter(1,7.876902);  
  fx->SetParameter(2,48.000000);   
  fx->SetParameter(3,44.100000);  
  fx->SetParameter(4,317.000000);  
  fx->SetParameter(5,0.111901);     
  fx->SetParameter(6,0.230769);   
  fx->SetParameter(7,204230.959837);
  fx->SetParameter(8,7140.115729);     

  std::vector<double> p;
  p.clear();
  for(int i=0;i<npar;i++) p.push_back(fx->GetParameter(i));

  //double compensation = th1->GetMinimum()-fx->GetParameter(0);
  double offset = -5600;
  TF1 *f0 = new TF1("f0","[0]",0,binnum*rebin/10);       // constant bg
  f0->SetParameter(0,p[7]);

  TF1 *f1 = new TF1("f1",Draw_fit,0,binnum*rebin/10,4);               // total bg
  f1->SetParameter(0,0);                                 // return parent curve
  f1->SetParameter(1,p[7]);                        // offset
  f1->SetParameter(2,p[8]);                        // slope
  f1->SetParameter(3,p[9]);                        // quad

  TF1 *f2 = new TF1("f2",Draw_fit,0,60,5);  // 30Ne
  f2->SetParameter(0,1);                                 // return parent curve
  f2->SetParameter(1,p[7]);                        // constant bg
  f2->SetParameter(2,1);                                 // braching ratio
  f2->SetParameter(3,p[0]);                        // actvity
  f2->SetParameter(4,p[1]);                        // parent t1/2
  TF1 *f22 = new TF1("f22",Draw_fit,0,binnum*rebin/10,6); // 30Ne
  f22->SetParameter(0,-1);                                 // return parent curve
  f22->SetParameter(1,p[7]);                         // A_bg
  f22->SetParameter(2,p[8]);                         // t1/2_bg
  f22->SetParameter(3,1);                                 // braching ratio
  f22->SetParameter(4,p[0]);                        // actvity
  f22->SetParameter(5,p[1]);                        // parent t1/2

  TF1 *f3 = new TF1("f3",Draw_fit,0,binnum*rebin/10,6);  // 30Na
  f3->SetParameter(0,2);                                 // model choose 
  f3->SetParameter(1,p[7]+offset);                        // constant bg
  f3->SetParameter(2,1-p[5]);                      // branching ratio
  f3->SetParameter(3,p[0]);                        // acitivity
  f3->SetParameter(4,p[1]);                        // parent(30Ne) t1/2
  f3->SetParameter(5,p[2]);                        // daughter(30Na) t1/2

  TF1 *f4 = new TF1("f4",Draw_fit,0,binnum*rebin/10,6);  // 29Na
  f4->SetParameter(0,2);                                 // return grand curve
  f4->SetParameter(1,p[7]+offset);                        // constant bg
  f4->SetParameter(2,p[5]);                        // braching ratio
  f4->SetParameter(3,p[0]);                        // activity
  f4->SetParameter(4,p[1]);                        // parent(30Ne) t1/2
  f4->SetParameter(5,p[3]);                        // daughter(29Na) t1/2

  TF1 *f5 = new TF1("f5",Draw_fit,0,binnum*rebin/10,7);  // 30Mg
  f5->SetParameter(0,3);                                 // return daughter curve
  f5->SetParameter(1,p[7]+offset);                        // constant bg
  f5->SetParameter(2,(1-p[5])*(1-p[6]));     // branching ratio
  f5->SetParameter(3,p[0]);                        // activity
  f5->SetParameter(4,p[1]);                        // parent(30Ne) t1/2
  f5->SetParameter(5,p[2]);                        // daughter(30Na) 1/2
  f5->SetParameter(6,p[4]);                        // grand(30Mg) t1/2

  

  f0->SetLineWidth(2);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  f22->SetLineColor(kGreen);
  f22->SetLineWidth(2);
  f3->SetLineColor(kViolet);
  f3->SetLineWidth(2);
  f4->SetLineColor(kCyan);
  f4->SetLineWidth(2);
  f5->SetLineColor(kBlack);
  f5->SetLineWidth(2);

  //gr->GetListOfFunctions()->Add(f0);
  gr->GetListOfFunctions()->Add(f1);
  gr->GetListOfFunctions()->Add(f22);
  gr->GetListOfFunctions()->Add(f3);
  gr->GetListOfFunctions()->Add(f4);
  gr->GetListOfFunctions()->Add(f5);

  //====================== convert TF1 to txt file ========================//
  std::ofstream ofile2;  
  ofile2.open("fit.dat");
  std::ofstream ofile3;  
  ofile3.open("exp_bg.dat");
  std::ofstream ofile4;  
  ofile4.open("fit30Ne.dat");
  std::ofstream ofile5;  
  ofile5.open("fit30Na.dat");
  std::ofstream ofile6;  
  ofile6.open("fit29Na.dat");
  std::ofstream ofile7;  
  ofile7.open("fit30Mg.dat");
  for(double i=0;i<1000;i=i+0.1){
    ofile2 << i << "\t" << fx->Eval(i)/100000. << std::endl;
    ofile3 << i << "\t" << f1->Eval(i)/100000. << std::endl;
    ofile4 << i << "\t" << f22->Eval(i)/100000. << std::endl;
    ofile5 << i << "\t" << f3->Eval(i)/100000. << std::endl;
    ofile6 << i << "\t" << f4->Eval(i)/100000. << std::endl;
    ofile7 << i << "\t" << f5->Eval(i)/100000. << std::endl;
  }

  return gr;

}




TGraph *simfithl(TH1D *hist1, int rebin=1, int binnum=-1){


  TH1D *th1 = (TH1D *)hist1->Clone("th1");
  th1->Rebin(rebin);
  double Ymax = th1->GetMaximum();
  int nbins = th1->GetNbinsX();
  if(binnum<0) binnum = nbins;
  if(binnum>nbins) { printf("WARNING: Fitting range is beyond the histogram!!!\n"); }

  int npar = 3;
  TGraph *gr = convert_th1_tgr(th1, binnum); 
  TF1 *fx = new TF1("fx","[0]*exp(-0.693/[1]*x)+[2]");

  fx->SetParameter(0,350);  // exp ion scaling factor
  fx->SetParameter(1,11.);   // ion half-life
  fx->SetParameter(2,0.1);   // constant bg                                               

  fx->FixParameter(0,225.629);
  fx->FixParameter(1,8.46414);
  fx->FixParameter(2,0.59972);
  //fx->SetParLimits(1,0,Ymax);    // exp ion scaling factor
  //fx->SetParLimits(2,0,10000);   // exp ion half-life

  gr->Fit(fx); 

  TF1 *f0 = new TF1("f0","[0]",0,100);     // constant bg
  f0->SetParameter(0,fx->GetParameter(2));

  TF1 *f1 = new TF1("f1","[0]*exp(-0.693/[1]*x)",0,100);    // 32Na
  f1->SetParameter(0,fx->GetParameter(0));    // actvity
  //f1->SetParameter(0,fx->GetParameter(1)+fx->GetParameter(0));    // actvity
  f1->SetParameter(1,fx->GetParameter(1));    // t1/2


  f0->SetLineWidth(2);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);

  gr->GetListOfFunctions()->Add(f0);
  gr->GetListOfFunctions()->Add(f1);

  //====================== convert TF1 to txt file ========================//
  std::ofstream ofile2;  
  ofile2.open("fit.dat");
  std::ofstream ofile3;  
  ofile3.open("flat_bg.dat");
  std::ofstream ofile4;  
  ofile4.open("exp_ion.dat");
  for(double i=0;i<150;i=i+0.1){
    ofile2 << i << "\t" << fx->Eval(i) << std::endl;
    ofile3 << i << "\t" << f0->Eval(i) << std::endl;
    ofile4 << i << "\t" << f1->Eval(i) << std::endl;
  }

  return gr;

}
