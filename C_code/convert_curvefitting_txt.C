

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

double SimFit(double *x, double *par){
  double value = par[0] + MYExp(x,par+1);
  return value;  
}

// Calculate activity of one decay = lambda *N(t);
// N(t) is calculated by Bateman equation;
// n = isotope number(parent = 1, daughter = 2, granddaughter = 3);
// when n = 1, return exponential function par[0]*exp(-lam1*x[0])
// par[0] = parent activity;par[i] = half-life of particle_i(i>0)
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
// n =  # of isotopes in this decay chain;
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
  //par[0]: constant bg;
  //par[1]: A_bg
  //par[2]: t1/2_bg;
  //par[3]: A_path1_32Na;
  //par[4]: t1/2_32Na; 
  //par[5]: t1/2_32Mg; 
  //par[6]: t1/2_32Al; 
  //par[7]: A_path2_32Na;
  //par[8]: t1/2_31Mg; 

  double par_bg[2];
  double par_path1[4];
  double par_path2[3];

  par_bg[0] = par[1];
  par_bg[1] = par[2];

  long nisotope1 = *(&par_path1+1)-par_path1 - 1; // array length -1(1st par is always activity)
  par_path1[0] = par[3];
  par_path1[1] = par[4];
  par_path1[2] = par[5];
  par_path1[3] = par[6];

  long nisotope2 = *(&par_path2+1)-par_path2 - 1;
  par_path2[0] = par[7];
  par_path2[1] = par[4];
  par_path2[2] = par[8];


  double value =  par[0] + MYExp(x,par_bg) + sum_Bateman(x,par_path1,3) + sum_Bateman(x,par_path2, 2);

  return value;
}
//================================== Math model for fitting plot =================================//
//par[0] = 0 => bg plot;
//par[0] >=0 => the value of par[0] = isotope number (parent=1, daughter=2, granddaughter=3)
double Draw_fit(double *x, double *par){
  double sum = 0;
  if(par[0]==0) sum = par[1] + MYExp(x,par+2);
  else sum = par[1] + Bateman(x,par+2,par[0]);
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
  for(int i=1;i<=binnum;i++){
    myx.push_back(th1->GetBinCenter(i));
    myy.push_back(th1->GetBinContent(i));
    mydx.push_back(0);                    // X-error
    mydy.push_back(sqrt(myy.back()));     // Y-error
    ofile1<<myx.back()<<"\t"
      <<myy.back()<<"\t"
      <<mydx.back()<<"\t"
      <<mydy.back()<<std::endl;
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
  TF1 *fx = new TF1("fx",MYFit,0,binnum,npar);

  fx->SetParameter(0,17500.4);   // constant bg                                               
  fx->SetParameter(1,9959.21);   // exp bg scaling factor
  fx->SetParameter(2,3751.9);   // exp bg half-life
  fx->SetParameter(3,7156.49);   // 32Na scaling factor in beta- path = 32Mg scaling factor
  fx->SetParameter(4,13.5323);  // 32Na half-life
  fx->SetParameter(5,86);     // 32Mg half-life
  fx->SetParameter(6,31.7);   // 32Al half-life
  fx->SetParameter(7,4768.62);   // 32Na scaling factor in beta-n path = 31Mg scaling factor
  fx->SetParameter(8,236);    // 31Mg half-life 

  fx->SetParLimits(0,0,Ymax);    // constant bg                                               
  fx->SetParLimits(1,0,Ymax);    // exp bg scaling factor
  fx->SetParLimits(2,0,10000);   // exp bg half-life
  fx->SetParLimits(3,0,Ymax);    // 32Na scaling factor in beta- path = 32Mg scaling factor
  fx->SetParLimits(4,10,16);     // 32Na half-life
  fx->SetParLimits(5,75,95);     // 32Mg half-life
  fx->SetParLimits(6,20,45);     // 32Al half-life
  fx->SetParLimits(7,0,Ymax);    // 32Na scaling factor in beta-n path = 31Mg scaling factor
  fx->SetParLimits(8,229,241);   // 31Mg half-life 

  fx->FixParameter(5,86.);
  fx->FixParameter(6,31.7);
  fx->FixParameter(8,236.);
  fx->SetLineWidth(2);
  gr->Fit(fx); 

  //double compensation = th1->GetMinimum()-fx->GetParameter(0);
  double compensation = 25750.-fx->GetParameter(0);

  TF1 *f0 = new TF1("f0","[0]",0,1000);     // constant bg
  f0->SetParameter(0,fx->GetParameter(0)+compensation);

  TF1 *f1 = new TF1("f1",Draw_fit,0,1000,4);  // 32Na
  f1->SetParameter(0,0);                      // return parent curve
  double A_tot = fx->GetParameter(3) + fx->GetParameter(7);
  f1->SetParameter(1,fx->GetParameter(0)+compensation);    // constant bg
  f1->SetParameter(2,A_tot);                // actvity
  f1->SetParameter(3,fx->GetParameter(4));    // t1/2

  TF1 *f2 = new TF1("f2",Draw_fit,0,1000,5);  // 32Mg
  f2->SetParameter(0,2);                      // return daughter curve
  f2->SetParameter(1,fx->GetParameter(0)+compensation);    // constant bg
  f2->SetParameter(2,fx->GetParameter(3));    // actvity
  f2->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  f2->SetParameter(4,fx->GetParameter(5));    // daughter t1/2

  TF1 *f3 = new TF1("f3",Draw_fit,0,1000,4);  // exp bg
  f3->SetParameter(0,0);                      // model choose 
  f3->SetParameter(1,fx->GetParameter(0));    // constant bg
  f3->SetParameter(2,fx->GetParameter(1));    // scaling factor
  f3->SetParameter(3,fx->GetParameter(2));    // t1/2

  TF1 *f4 = new TF1("f4",Draw_fit,0,1000,6);  // 32Al
  f4->SetParameter(0,3);                      // return grand curve
  f4->SetParameter(1,fx->GetParameter(0)+compensation);    // constant bg
  f4->SetParameter(2,fx->GetParameter(3));    // activity
  f4->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  f4->SetParameter(4,fx->GetParameter(5));    // daughter t1/2
  f4->SetParameter(5,fx->GetParameter(6));    // grand t1/2

  TF1 *f5 = new TF1("f5",Draw_fit,0,1000,5);  // 31Mg
  f5->SetParameter(0,2);                      // return daughter curve
  f5->SetParameter(1,fx->GetParameter(0)+compensation);    // constant bg
  f5->SetParameter(2,fx->GetParameter(7));    // activity_2
  f5->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  f5->SetParameter(4,fx->GetParameter(8));    // daughter 1/2

  f0->SetLineWidth(2);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(2);
  f3->SetLineColor(kViolet);
  f3->SetLineWidth(2);
  f4->SetLineColor(kCyan);
  f4->SetLineWidth(2);
  f5->SetLineColor(kBlack);
  f5->SetLineWidth(2);

  gr->GetListOfFunctions()->Add(f0);
  gr->GetListOfFunctions()->Add(f1);
  gr->GetListOfFunctions()->Add(f2);
  gr->GetListOfFunctions()->Add(f3);
  gr->GetListOfFunctions()->Add(f4);
  gr->GetListOfFunctions()->Add(f5);

  //====================== convert TF1 to txt file ========================//
  std::ofstream ofile2;  
  ofile2.open("fit.dat");
  std::ofstream ofile3;  
  ofile3.open("flat_bg.dat");
  std::ofstream ofile4;  
  ofile4.open("exp_bg.dat");
  std::ofstream ofile5;  
  ofile5.open("fit32Na.dat");
  std::ofstream ofile6;  
  ofile6.open("fit32Mg.dat");
  std::ofstream ofile7;  
  ofile7.open("fit32Al.dat");
  std::ofstream ofile8;  
  ofile8.open("fit31Mg.dat");
  for(double i=0;i<1000;i=i+0.1){
    ofile2 << i << "\t" << fx->Eval(i) << std::endl;
    ofile3 << i << "\t" << f0->Eval(i) << std::endl;
    ofile4 << i << "\t" << f3->Eval(i) << std::endl;
    ofile5 << i << "\t" << f1->Eval(i) << std::endl;
    ofile6 << i << "\t" << f2->Eval(i) << std::endl;
    ofile7 << i << "\t" << f4->Eval(i) << std::endl;
    ofile8 << i << "\t" << f5->Eval(i) << std::endl;
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
