

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


/*double Bateman(double *x, double *par, int n=1){
  
  double val = 0;
  double lam1, lam2, lam3;
  if(n==1){
    lam1 = 0.693/par[1];
    val = par[0]*TMath::Exp(-lam1*x[0]);
  }
  if(n==2){
    lam1 = 0.693/par[1];
    lam2 = 0.693/par[2];
    val = lam2*par[0]/(lam2-lam1)*(TMath::Exp(-lam1*x[0])-TMath::Exp(-lam2*x[0]));
  }
  if(n==3){
    lam1 = 0.693/par[1];
    lam2 = 0.693/par[2];
    lam3 = 0.693/par[3];
    val = lam3*par[0]*lam2*(TMath::Exp(-0.693/par[1]*x[0])/(lam2-lam1)/(lam3-lam1)+TMath::Exp(-0.693/par[2]*x[0])/(lam1-lam2)/(lam3-lam2)+TMath::Exp(-0.693/par[3]*x[0])/(lam1-lam3)/(lam2-lam3));
  }
  return val;

}*/

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
  //par[0]: constant bg;
  //par[1]: A_bg
  //par[2]: t1/2_bg;
  //par[3]: A: path1: 30Ne->30Na->30Mg->30Al;
  //par[4]: A: path2: 30Ne->30Na->29Mg;
  //par[5]: A: path3: 30Ne->29Na->29Mg
  //par[6]: t1/2: 30Ne;
  //par[7]: t1/2: 30Na;   
  //par[8]: t1/2: 29Na;
  //par[9]: t1/2: 30Mg;

 
  double par_bg[2];
  double par1[3]; //path1: 30Ne->30Na->30Mg->30Al
  double par2[2]; //path2: 30Ne->30Na->29Mg
  double par3[3]; //path3: 30Ne->29Na->29Mg

  par_bg[0] = par[1];
  par_bg[1] = par[2];

  par1[0] = par[3];
  par1[1] = par[6];
  par1[2] = par[7];
  //par1[3] = par[9];

  par2[0] = par[4];
  par2[1] = par[6];
  //par2[2] = par[7];
  
  par3[0] = par[5];
  par3[1] = par[6];
  par3[2] = par[8];
  
  double value = par[0] + MYExp(x,par_bg) + sum_Bateman(x,par1,2);// + sum_Bateman(x,par2,2);// + sum_Bateman(x,par3,2);
  
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
std::pair<std::vector<double>, std::vector<double>> searchinitPar(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2.){
  
  int npar = 10;
  std::pair<std::vector<double>,std::vector<double>> presult;  
  
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
  p[1] = 1000/20.*rebin;    //A_bg
  p[2] = 340;               //t1/2_bg;
  p[3] = 5000/20.*rebin;    //A: 30Ne->30Na->30Mg->30Al
  p[4] = 500/20.*rebin;     //A: 30Ne->30Na->29Mg;
  p[5] = 500/20.*rebin;     //A: 30Ne->29Na->29Mg;
  p[6] = 7.3;               //t1/2_30Ne;
  p[7] = 44;                //t1/2_30Na
  p[8] = 40;                //t1/2_29Na;
  p[9] = 300;               //t1/2_30Mg;

  TF1 *fx = new TF1("fx",MYFit,0,binnum*rebin/10,npar);
  for(int i=0;i<npar;i++){
    fx->SetParameter(i,p[i]);
  }
  fx->FixParameter(1,0);
  fx->FixParameter(2,0);
  //fx->FixParameter(3,0);
  fx->FixParameter(4,0);
  fx->FixParameter(5,0);
  //fx->FixParameter(7,0);
  fx->FixParameter(8,0);
  fx->FixParameter(9,0);

  gr->Fit(fx); 
  double chi2 = fx->GetChisquare();
  double ndf = fx->GetNDF();
  for(int i=0;i<npar;i++){
    presult.first.push_back(fx->GetParameter(i));
  }

  // change parmeter guess linealy(You can change math formular and conditions to stop)
  int iteration=0;
  while(chi2>1){
    for(int ipar=0;ipar<npar;ipar++){
      fx->SetParameter(ipar,presult.first[ipar]);
    }
    fx->SetParLimits(0,Ymin,Ymax);     //constant bg                                          
    fx->SetParLimits(1,0,Ymax);        //A_bg
    fx->SetParLimits(2,0,10000);       //t1/2_bg;
    fx->SetParLimits(3,1000,Ymax);     //A: 30Ne->30Na->30Mg->30Al
    fx->SetParLimits(4,1,Ymax);        //A: 30Ne->30Na->29Mg
    fx->SetParLimits(5,1,Ymax);        //A: 30Ne->29Na->29Mg
    fx->SetParLimits(6,5,9);       //t1/2: 30Ne;
    fx->SetParLimits(7,40,52);         //t1/2: 30Na
    fx->SetParLimits(8,42,44);         //t1/2: 29Na
    fx->SetParLimits(9,60,330);       //t1/2: 30Mg;
    
    fx->FixParameter(1,0);
    fx->FixParameter(2,0);
    //fx->FixParameter(3,0);
    fx->FixParameter(4,0);
    fx->FixParameter(5,0);
    //fx->FixParameter(7,0);
    fx->FixParameter(8,0);
    fx->FixParameter(9,317);

    gr->Fit(fx);

    double test_chi2 = fx->GetChisquare();
    test_chi2 = test_chi2/ndf;
    // if chi2<min_chi2 && chi2 start increasing OR the current chi2<1 OR iteration over 100times,
    // stop searching initial parameters and return;
    if((test_chi2 > chi2 && test_chi2<chi2_min) || (test_chi2<1) || iteration>1000){
      if(fabs(chi2-1)>fabs(test_chi2-1) && chi2>2){
        printf("chi2 = %f\n", test_chi2);
        presult.first.clear();
        presult.second.clear();
        for(int i=0;i<npar;i++){
          presult.first.push_back(fx->GetParameter(i));
          presult.second.push_back(fx->GetParError(i));
          printf("p[%i] = %f\n", i,presult.first[i]);
        }
      }else{
        printf("chi2 = %f\n", chi2);
        for(int i=0;i<npar;i++){
          printf("p[%i] = %f\n", i,presult.first[i]);
        }
      }
      break;
    }
    if(test_chi2<chi2){
      presult.first.clear();
      presult.second.clear();
      for(int i=0;i<npar;i++){
        presult.first.push_back(fx->GetParameter(i));
        presult.second.push_back(fx->GetParError(i));
      }
      presult.first.push_back(test_chi2);
    }
    chi2 = test_chi2;
    iteration++;
  }

  return presult;

}

//============================= Draw fitting results ============================//

void FitHL(TH1D *hist1, int rebin=1, int binnum=-1, double chi2_min=2){

  std::pair<std::vector<double>, std::vector<double>> p = searchinitPar(hist1, rebin, binnum, chi2_min);
  int npar = p.first.size()-1;
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
    fx->SetParameter(i,p.first[i]);
  }
  fx->SetLineColor(kRed);
  fx->SetLineWidth(2);
  gr->GetListOfFunctions()->Add(fx);
  
  TF1 *f0 = new TF1("f0","[0]",0,binnum*rebin/10);     // constant bg
  f0->SetParameter(0,p.first[0]);

  TF1 *f1 = new TF1("f1",Draw_fit,0,binnum*rebin/10,4);  // exp bg
  f1->SetParameter(0,1);                                 // return parent curve
  f1->SetParameter(1,p.first[0]);                        // constant bg
  f1->SetParameter(2,p.first[1]);                        // actvity
  f1->SetParameter(3,p.first[2]);                        // t1/2

  TF1 *f2 = new TF1("f2",Draw_fit,0,binnum*rebin/10,4);  // 30Ne
  f2->SetParameter(0,1);                      // return parent curve
  f2->SetParameter(1,p.first[0]);             // constant bg
  f2->SetParameter(2,p.first[3]+p.first[4]+p.first[5]);             // actvity
  f2->SetParameter(3,p.first[6]);             // parent t1/2

  TF1 *f3 = new TF1("f3",Draw_fit,0,binnum*rebin/10,5);  // 30Na
  f3->SetParameter(0,2);                      // model choose 
  f3->SetParameter(1,p.first[0]);             // constant bg
  f3->SetParameter(2,p.first[3]+p.first[4]);             // acitivity
  f3->SetParameter(3,p.first[6]);             // parent(30Ne) t1/2
  f3->SetParameter(4,p.first[7]);             // daughter(30Na) t1/2

  TF1 *f4 = new TF1("f4",Draw_fit,0,binnum*rebin/10,5);  // 29Na
  f4->SetParameter(0,2);                      // return grand curve
  f4->SetParameter(1,p.first[0]);             // constant bg
  f4->SetParameter(2,p.first[5]);             // activity
  f4->SetParameter(3,p.first[6]);             // parent(30Ne) t1/2
  f4->SetParameter(4,p.first[8]);             // daughter(29Na) t1/2

  TF1 *f5 = new TF1("f5",Draw_fit,0,binnum*rebin/10,6);  // 30Mg
  f5->SetParameter(0,3);                      // return daughter curve
  f5->SetParameter(1,p.first[0]);             // constant bg
  f5->SetParameter(2,p.first[3]);             // activity
  f5->SetParameter(3,p.first[6]);             // parent(30Ne) t1/2
  f5->SetParameter(4,p.first[7]);             // daughter(30Na) 1/2
  f5->SetParameter(5,p.first[9]);             // grand(30Mg) t1/2

  f0->SetLineWidth(2);
  f0->SetTitle("flat bg");
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  f1->SetTitle("exp bg");
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(2);
  f2->SetTitle("30Ne");
  f3->SetLineColor(kViolet);
  f3->SetLineWidth(2);
  f3->SetTitle("30Na");
  f4->SetLineColor(kCyan);
  f4->SetLineWidth(2);
  f4->SetTitle("29Na");
  f5->SetLineColor(kBlack);
  f5->SetLineWidth(2);
  f5->SetTitle("30Mg");

  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  gr->GetListOfFunctions()->Add(f0);
  gr->GetListOfFunctions()->Add(f1);
  gr->GetListOfFunctions()->Add(f2);
  gr->GetListOfFunctions()->Add(f3);
  gr->GetListOfFunctions()->Add(f4);
  gr->GetListOfFunctions()->Add(f5);

  gr->GetYaxis()->SetRangeUser(p.first[0]*0.995, Ymax*1.05);

  TCanvas *c = new TCanvas;
  gr->Draw();
  
  auto legend = new TLegend(0.2,0.7,0.4,0.9);
  legend->SetNColumns(2);
  legend->AddEntry(f0, "flat bg", "l");
  legend->AddEntry(f1, "exp bg", "l");
  legend->AddEntry(f2, "30Ne", "l");
  legend->AddEntry(f3, "30Na", "l");
  legend->AddEntry(f4, "29Na", "l");
  legend->AddEntry(f5, "30Mg", "l");
  legend->Draw();
  

}

void hlerrcheck(TH1D *hist, int rebin, double chi2_min=-1){
 
  
  int size = 9;
  int binnum;
  if(chi2_min<0) chi2_min=2;
  TGraph *tempgr[size];
  std::pair<std::vector<double>, std::vector<double>> p[size];
  TH1D *th1[size]; 
  TF1 *fx[size];
  TF1 *f0[size];
  TF1 *f1[size];
  TF1 *f2[size];
  TF1 *f3[size];
  TF1 *f4[size];
  TF1 *f5[size];
  TLegend *legend[size];
   
  TCanvas *c = new TCanvas;
  c->SetTitle(Form("bin in %.1fms",rebin/10.));
  c->Divide(3,3);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  for(int i=0;i<size;i++){
    binnum = (1500+500*i)/rebin;
    p[i].first.clear();
    p[i].second.clear();
    th1[i] = (TH1D *)hist->Clone(Form("decaytime_%ims_rebin_%i",binnum/10,rebin));
    th1[i]->Rebin(rebin);
    double Ymax = th1[i]->GetMaximum();
    double Ymin = th1[i]->GetMinimum();
    tempgr[i] = convert_th1_tgr(th1[i], binnum);
    tempgr[i]->SetTitle(Form("decaytime = %ims, bin in %.1fms",binnum*rebin/10,rebin/10.));
    p[i] = searchinitPar(hist,rebin,binnum,chi2_min);
    fx[i] = new TF1(Form("fx%i",i),MYFit,0,binnum*rebin/10,p[i].first.size()-1);
    fx[i]->SetLineColor(kRed);
    fx[i]->SetLineWidth(2);
    for(int npar=0;npar<p[i].second.size();npar++){
      fx[i]->SetParameter(npar,p[i].first[npar]);
    } 
    tempgr[i]->GetListOfFunctions()->Add(fx[i]);
    f0[i] = new TF1(Form("f0_%i",i),"[0]",0,binnum*rebin/10);     // constant bg
    f0[i]->SetParameter(0,fx[i]->GetParameter(0));

    f1[i] = new TF1(Form("f1_%i",i),Draw_fit,0,binnum*rebin/10,4);  // exp bg
    f1[i]->SetParameter(0,1);                                 // return parent curve
    f1[i]->SetParameter(1,p[i].first[0]);                        // constant bg
    f1[i]->SetParameter(2,p[i].first[1]);                        // actvity
    f1[i]->SetParameter(3,p[i].first[2]);                        // t1/2

    f2[i] = new TF1(Form("f2_%i",i),Draw_fit,0,binnum*rebin/10,4);  // 30Ne
    f2[i]->SetParameter(0,1);                      // return parent curve
    f2[i]->SetParameter(1,p[i].first[0]);             // constant bg
    f2[i]->SetParameter(2,p[i].first[3]+p[i].first[4]+p[i].first[5]);             // actvity
    f2[i]->SetParameter(3,p[i].first[6]);             // parent t1/2

    f3[i] = new TF1(Form("f3_%i",i),Draw_fit,0,binnum*rebin/10,5);  // 30Na
    f3[i]->SetParameter(0,2);                      // model choose 
    f3[i]->SetParameter(1,p[i].first[0]);             // constant bg
    f3[i]->SetParameter(2,p[i].first[3]+p[i].first[4]);             // acitivity
    f3[i]->SetParameter(3,p[i].first[6]);             // parent(30Ne) t1/2
    f3[i]->SetParameter(4,p[i].first[7]);             // daughter(30Na) t1/2

    f4[i] = new TF1(Form("f4_%i",i),Draw_fit,0,binnum*rebin/10,5);  // 29Na
    f4[i]->SetParameter(0,2);                      // return grand curve
    f4[i]->SetParameter(1,p[i].first[0]);             // constant bg
    f4[i]->SetParameter(2,p[i].first[5]);             // activity
    f4[i]->SetParameter(3,p[i].first[6]);             // parent(30Ne) t1/2
    f4[i]->SetParameter(4,p[i].first[8]);             // daughter(29Na) t1/2

    f5[i] = new TF1(Form("f5_%i",i),Draw_fit,0,binnum*rebin/10,6);  // 30Mg
    f5[i]->SetParameter(0,3);                      // return daughter curve
    f5[i]->SetParameter(1,p[i].first[0]);             // constant bg
    f5[i]->SetParameter(2,p[i].first[3]);             // activity
    f5[i]->SetParameter(3,p[i].first[6]);             // parent(30Ne) t1/2
    f5[i]->SetParameter(4,p[i].first[7]);             // daughter(29Na) 1/2
    f5[i]->SetParameter(5,p[i].first[9]);             // grand(30Mg) t1/2

    f0[i]->SetLineWidth(2);
    f0[i]->SetTitle("flat bg");
    f1[i]->SetLineColor(kBlue);
    f1[i]->SetLineWidth(2);
    f1[i]->SetTitle("exp bg");
    f2[i]->SetLineColor(kGreen);
    f2[i]->SetLineWidth(2);
    f2[i]->SetTitle("30Ne");
    f3[i]->SetLineColor(kViolet);
    f3[i]->SetLineWidth(2);
    f3[i]->SetTitle("30Na");
    f4[i]->SetLineColor(kCyan);
    f4[i]->SetLineWidth(2);
    f4[i]->SetTitle("29Na");
    f5[i]->SetLineColor(kBlack);
    f5[i]->SetLineWidth(2);
    f5[i]->SetTitle("30Mg");

    tempgr[i]->GetListOfFunctions()->Add(f0[i]);
    tempgr[i]->GetListOfFunctions()->Add(f1[i]);
    tempgr[i]->GetListOfFunctions()->Add(f2[i]);
    tempgr[i]->GetListOfFunctions()->Add(f3[i]);
    tempgr[i]->GetListOfFunctions()->Add(f4[i]);
    tempgr[i]->GetListOfFunctions()->Add(f5[i]);

    tempgr[i]->GetYaxis()->SetRangeUser(p[i].first[0]*0.995, Ymax*1.003);
    

    c->cd(i+1);
    tempgr[i]->Draw();
 
    legend[i] = new TLegend(0.7,0.7,0.9,0.9);
    //legend[i] = new TLegend();
    legend[i]->AddEntry(f0[i], "flat bg", "l");
    legend[i]->AddEntry(f1[i], "exp bg", "l");
    legend[i]->AddEntry(f2[i], "30Ne", "l");
    legend[i]->AddEntry(f3[i], "30Na", "l");
    legend[i]->AddEntry(f4[i], "29Na", "l");
    legend[i]->AddEntry(f5[i], "30Mg", "l");
    //legend[i]->Draw();
 
    double y = gPad->GetUymax();
    TText ptext[11];
    for(int j=0;j<p[i].first.size();j++){
      ptext[j].SetNDC();
      ptext[j].SetTextFont(1);
      ptext[j].SetTextColor(1);
      ptext[j].SetTextSize(0.05);
      if(j<p[i].first.size()-1){
        ptext[j].DrawText(0.3, y*(0.85-0.05*j), Form("p[%i] = %f(%f)",j,p[i].first[j],p[i].second[j]));
      }else{
        ptext[j].DrawText(0.3, y*(0.85-0.05*j), Form("chi2 = %f",p[i].first[j]));
      }
    }
    //TText *chi2 = new TText();
    //chi2 -> SetNDC();
    //chi2 -> SetTextFont(1);
    //chi2 -> SetTextColor(1);
    //chi2 -> SetTextSize(0.05);
    //chi2 -> DrawText(0.3, y*.8, Form("X2 = %.2f",p[i].first.back()));
    //TText *hl = new TText();
    //hl -> SetNDC();
    //hl -> SetTextFont(1);
    //hl -> SetTextColor(1);
    //hl -> SetTextSize(0.05);
    //if(p[i].second.size()<5){
    //  hl -> DrawText(0.3, y*.75, Form("hl = %.2fms",p[i].first[6]));
    //}else{
    //  hl -> DrawText(0.3, y*.75, Form("hl = %.2f(%.2f)ms",p[i].first[6], p[i].second[6]));
    //} 
    //TText *Act = new TText();
    //Act -> SetNDC();
    //Act -> SetTextFont(1);
    //Act -> SetTextColor(1);
    //Act -> SetTextSize(0.05);
    //Act -> DrawText(0.3, y*.7, Form("A(30Ne) = %f",p[i].first[3]+p[i].first[4]+p[i].first[5]));
    //TText *Act1 = new TText();
    //Act1 -> SetNDC();
    //Act1 -> SetTextFont(1);
    //Act1 -> SetTextColor(1);
    //Act1 -> SetTextSize(0.05);
    //Act1 -> DrawText(0.3, y*.65, Form("A(30Ne->30Na) = %f",p[i].first[3]+p[i].first[4]));
    //TText *Act2 = new TText();
    //Act2 -> SetNDC();
    //Act2 -> SetTextFont(1);
    //Act2 -> SetTextColor(1);
    //Act2 -> SetTextSize(0.05);
    //Act2 -> DrawText(0.3, y*.6, Form("A(30Ne->29Na) = %f",p[i].first[4]));
    //TText *Act3 = new TText();
    //Act3 -> SetNDC();
    //Act3 -> SetTextFont(1);
    //Act3 -> SetTextColor(1);
    //Act3 -> SetTextSize(0.05);
    //Act3 -> DrawText(0.3, y*.55, Form("A(30Ne->30Mg) = %f",p[i].first[3]));
    //TText *hlbg = new TText();
    //hlbg -> SetNDC();
    //hlbg -> SetTextFont(1);
    //hlbg -> SetTextColor(1);
    //hlbg -> SetTextSize(0.05);
    //hlbg -> DrawText(0.3, y*.6, Form("t1/2(bg) = %.2f",p[i].first[2]));
    

    c->Update();

  }

}


