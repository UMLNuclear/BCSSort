

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
  //par[3]: A_path1_30Ne
  //par[4]: t1/2_30Ne;
  //par[5]: t1/2_30Na
  //par[6]: A_path2_30Ne;
  //par[7]: A_path3_30Ne;
  //par[8]: t1/2_29Na
  //par[9]: t1/2_30Mg
  
  double par_bg[2];
  double par1[4]; //30Ne->30Na->30Mg->30Al
  double par2[3]; //30Ne->30Na->29Mg
  double par3[3]; //30Ne->29Na->29Mg

  par_bg[0] = par[1];
  par_bg[1] = par[2];

  par1[0] = par[3];
  par1[1] = par[4];
  par1[2] = par[5];
  par1[3] = par[9];

  par2[0] = par[6];
  par2[1] = par[4];
  par2[2] = par[5];

  par3[0] = par[7];
  par3[1] = par[4];
  par3[2] = par[8];

  double value = par[0] + MYExp(x,par_bg) + sum_Bateman(x,par1,3) + sum_Bateman(x,par2,2) + sum_Bateman(x,par3,2);
  
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

  TGraphErrors *gr = 0;
  std::vector<double> myx;
  std::vector<double> myy;
  std::vector<double> mydx;
  std::vector<double> mydy;
  for(int i=2;i<=binnum;i++){ //i=2 to skip 1st bin
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

//===================== return fitting parameters with min chi2 =================================//
std::vector<double> searchinitPar(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2.){
//TGraph* searchinitPar(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2.){
  int npar = 10;
  std::vector<double> presult;  
  
  TH1D *th1 = (TH1D *)hist1->Clone("th1");
  th1->Rebin(rebin);
  int nbins = th1->GetNbinsX();
  if(binnum<0) binnum = nbins;
  if(binnum>nbins) { printf("WARNING: Fitting range is beyond the histogram!!!\n"); }
  TGraph *gr = convert_th1_tgr(th1, binnum); 
  double Ymax = th1->GetMaximum();
  double Ymin = th1->GetMinimum();

  // Initila Pars Guess
  double p[npar];
  p[0] = 50000./20.*rebin;  //constant bg;
  p[1] = 1000/20.*rebin;   //A_bg
  p[2] = 340;               //t1/2_bg;
  p[3] = 5000/20.*rebin;    //A_path1_30Ne
  p[4] = 7.3;               //t1/2_30Ne;
  p[5] = 40;                //t1/2_30Na
  p[6] = 600/20.*rebin;     //A_path2_30Ne;
  p[7] = 250/20.*rebin;     //A_path3_30Ne;
  p[8] = 44.1;              //t1/2_29Na;
  p[8] = 317;               //t1/2_30Mg;

  TF1 *fx = new TF1("fx",MYFit,0,binnum,npar);
  for(int i=0;i<npar;i++){
    fx->SetParameter(i,p[i]);
  }
  fx->FixParameter(5,48);
  fx->FixParameter(8,44.1);
  fx->FixParameter(9,317);

  gr->Fit(fx,"N"); 
  double chi2 = fx->GetChisquare();
  double ndf = fx->GetNDF();
  for(int i=0;i<npar;i++){
    presult.push_back(fx->GetParameter(i));
  }

  // change parmeter guess linealy(You can change math formular and conditions to stop)
  while(chi2>1){
    p[0] = p[0] + 10;
    if(p[0]>Ymax) {printf("p[0] = %f\n",p[0]); break;}
    fx->SetParameter(0,p[0]);  
    fx->SetParLimits(0,Ymin,Ymax); 

    p[1] = p[1] + 10;
    if(p[1]>Ymax) {printf("p[0] = %f\n",p[1]); break;}
    fx->SetParameter(1,p[1]);  
    fx->SetParLimits(1,0,Ymax); 
   
    p[2] = p[2] + 10;
    if(p[2]>1500) {printf("p[2] = %f\n",p[2]); break;}
    fx->SetParameter(2,p[2]);
    fx->SetParLimits(2,0,10000);

    p[3] = p[3] +10;
    if(p[3]>Ymax) {printf("p[7] = %f\n",p[3]); break;}
    fx->SetParameter(3,p[3]);  
    fx->SetParLimits(3,1000,Ymax); 
       
    p[4] = p[4] + 0.05;
    if(p[4]>8.3){
      p[4] = 6.3;
      printf("reset p[4] = 6\n");
    }
    fx->SetParameter(4,p[4]);  
    fx->SetParLimits(4,6.3,8.3); 
    
    p[5] = p[5] + 0.25;
    if(p[5]>60){
      p[5] = 35;
      printf("reset p[5] = 35\n");
    }
    fx->SetParameter(5,p[5]);  
    fx->SetParLimits(5,35,55); 

    p[6] = p[6] +10;
    if(p[6]>Ymax) {printf("p[6] = %f\n",p[6]); break;}
    fx->SetParameter(6,p[6]);  
    fx->SetParLimits(6,50,Ymax); 

    p[7] = p[7] +10;
    if(p[7]>Ymax) {printf("p[7] = %f\n",p[7]); break;}
    fx->SetParameter(7,p[7]);  
    fx->SetParLimits(7,0,Ymax); 
 
    p[8] = p[8] + 0.25;
    if(p[8]>50){
      p[8] = 39;
      printf("reset p[8] = 39\n");
    }
    fx->SetParameter(8,p[8]);  
    fx->SetParLimits(8,40,50);
    
    p[9] = p[9] + 0.5;
    if(p[9]>330){
      p[9] = 310;
      printf("reset p[9] = 310\n");
    }
    fx->SetParameter(9,p[9]);  
    fx->SetParLimits(9,310,325);
    
    fx->FixParameter(5,48);
    fx->FixParameter(8,44.1);
    fx->FixParameter(9,317);

    gr->Fit(fx);

    double test_chi2 = fx->GetChisquare();
    test_chi2 = test_chi2/ndf;
    // if chi2<min_chi2 && chi2 start increasing OR the current chi2<1,
    // stop searching initial parameters and return;
    if((test_chi2 > chi2 && test_chi2<chi2_min) || (test_chi2<1)){
      printf("chi2 = %f\n", chi2);
      for(int i=0;i<npar;i++){
        printf("p[%i] = %f\n", i,presult[i]);
      }
      return presult;
      break;
    }
    chi2 = test_chi2;
    presult.clear();
    for(int i=0;i<npar;i++){
      presult.push_back(fx->GetParameter(i));
      //presult.push_back(p[i]);
    }
    presult.push_back(chi2);
     
  }

  return presult;

}

//============================= Draw fitting results ============================//

TGraph *FitHL(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2){

  std::vector<double> p = searchinitPar(hist1, rebin, binnum, chi2_min);
  int npar = p.size();
  for(int i=0;i<npar;i++){
    printf("FitHL p[%i] = %f\n", i,p[i]);
  }
  TH1D *th1 = (TH1D *)hist1->Clone("th1");
  th1->Rebin(rebin);
  double Ymax = th1->GetMaximum();
  double Ymin = th1->GetMinimum();
  int nbins = th1->GetNbinsX();
  if(binnum<0) binnum = nbins;
  if(binnum>nbins) { printf("WARNING: Fitting range is beyond the histogram!!!\n"); }

  TGraph *gr = convert_th1_tgr(th1, binnum); 
  TF1 *fx = new TF1("fx",MYFit,0,binnum*rebin/10,npar);
  for(int i=0;i<npar;i++){
    fx->SetParameter(i,p[i]);
  }
  fx->SetLineColor(kRed);
  fx->SetLineWidth(2);
  gr->GetListOfFunctions()->Add(fx);
  //fx->SetParLimits(0,Ymin,Ymax);     //constant bg                                          
  //fx->SetParLimits(1,0,Ymax);      //A_bg
  //fx->SetParLimits(2,0,10000);       //t1/2_bg;
  //fx->SetParLimits(3,1000,Ymax);        //A_path1_30Ne
  //fx->SetParLimits(4,6.3,8.3);       //t1/2_30Ne;
  //fx->SetParLimits(5,35,55);         //t1/2_30Na
  //fx->SetParLimits(6,50,Ymax);        //A_path2_30Ne
  //fx->SetParLimits(7,0,Ymax);        //A_path3_30Ne
  //fx->SetParLimits(8,40,50);         //t1/2_29Na;
  //fx->SetParLimits(9,310,325);       //t1/2_30Mg;

  //fx->FixParameter(5,48);
  //fx->FixParameter(8,44.1);
  //fx->FixParameter(9,317);

  //gr->Fit(fx); 

  //TF1 *f0 = new TF1("f0","[0]",0,1000);     // constant bg
  //f0->SetParameter(0,fx->GetParameter(0));

  //TF1 *f1 = new TF1("f1",Draw_fit,0,1000,4);  // 32Na
  //f1->SetParameter(0,1);                      // return parent curve
  //double A_tot = fx->GetParameter(3) + fx->GetParameter(7);
  //f1->SetParameter(1,fx->GetParameter(0));    // constant bg
  //f1->SetParameter(2,A_tot);                // actvity
  //f1->SetParameter(3,fx->GetParameter(4));    // t1/2

  //TF1 *f2 = new TF1("f2",Draw_fit,0,1000,5);  // 32Mg
  //f2->SetParameter(0,2);                      // return daughter curve
  //f2->SetParameter(1,fx->GetParameter(0));    // constant bg
  //f2->SetParameter(2,fx->GetParameter(3));    // actvity
  //f2->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  //f2->SetParameter(4,fx->GetParameter(5));    // daughter t1/2

  //TF1 *f3 = new TF1("f3",Draw_fit,0,1000,4);  // exp bg
  //f3->SetParameter(0,0);                      // model choose 
  //f3->SetParameter(1,fx->GetParameter(0));    // constant bg
  //f3->SetParameter(2,fx->GetParameter(1));    // scaling factor
  //f3->SetParameter(3,fx->GetParameter(2));    // t1/2

  //TF1 *f4 = new TF1("f4",Draw_fit,0,1000,6);  // 32Al
  //f4->SetParameter(0,3);                      // return grand curve
  //f4->SetParameter(1,fx->GetParameter(0));    // constant bg
  //f4->SetParameter(2,fx->GetParameter(3));    // activity
  //f4->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  //f4->SetParameter(4,fx->GetParameter(5));    // daughter t1/2
  //f4->SetParameter(5,fx->GetParameter(6));    // grand t1/2

  //TF1 *f5 = new TF1("f5",Draw_fit,0,1000,5);  // 31Mg
  //f5->SetParameter(0,2);                      // return daughter curve
  //f5->SetParameter(1,fx->GetParameter(0));    // constant bg
  //f5->SetParameter(2,fx->GetParameter(7));    // activity_2
  //f5->SetParameter(3,fx->GetParameter(4));    // parent t1/2
  //f5->SetParameter(4,fx->GetParameter(8));    // daughter 1/2

  //f0->SetLineWidth(2);
  //f1->SetLineColor(kBlue);
  //f1->SetLineWidth(2);
  //f2->SetLineColor(kGreen);
  //f2->SetLineWidth(2);
  //f3->SetLineColor(kViolet);
  //f3->SetLineWidth(2);
  //f4->SetLineColor(kCyan);
  //f4->SetLineWidth(2);
  //f5->SetLineColor(kBlack);
  //f5->SetLineWidth(2);

  //gr->GetListOfFunctions()->Add(f0);
  //gr->GetListOfFunctions()->Add(f1);
  //gr->GetListOfFunctions()->Add(f2);
  //gr->GetListOfFunctions()->Add(f3);
  //gr->GetListOfFunctions()->Add(f4);
  //gr->GetListOfFunctions()->Add(f5);


  return gr;

}


