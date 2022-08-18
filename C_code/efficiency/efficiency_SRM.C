#ifndef __CINT__
#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>


//#include <TApplication.h>
#include <TRint.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TString.h>

#endif




//==================== GammaEff ==============//
//Efficieny function
double GammaEff(double *x, double *par){
  //double logE = TMath::Log(x[0]);
  //double temp =  par[0] + par[1]*logE + par[2]*logE*logE + par[3]*pow(logE,3) +par[4]/(x[0]*x[0]);
  double logE = TMath::Log10(x[0]);
  double temp =  par[0] + par[1]*logE + par[2]*logE*logE + par[3]/(x[0]*x[0]);
  return pow(10,temp);
  //return pow(TMath::E(),temp);
}

//============================= NegRemove ==================//
//Remove bins with negative contents
//Return TH1 *
TH1 *NegRemove(TH1 *input){

  double nbins = input->GetNbinsX();
  double width = input->GetBinWidth(1);
  double lower = input->GetBinCenter(1);
  lower -= width/2;
  double upper = input->GetBinCenter(nbins);
  upper += width/2;
  TH1D *output = new TH1D(Form("%s_negremove",input->GetName()), Form("NegRmove %s",input->GetTitle()), nbins,lower,upper);
  for(int i=1;i<=nbins;i++){
    if(input->GetBinContent(i)<0){
      output->SetBinContent(i,0);
    }else{
      output->SetBinContent(i,input->GetBinContent(i));
    }
  }

  return output;
}


//========================================== SRMSinglePeakFit ================================//
//Fit one peak from SRM with different fitting range
//Draw fitting after bgsub with and without negative bin contents
//Input en = the number of energy:
//0 = 86keV;
//1 = 105keV;
//2 = 123keV;
void SRMSinglePeakFit(int en = 0, int xtal = 0, int iteration=10){
  
  double energy[10] = {  86.5, 105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[3] = {82, 103.5, 118.5};
  double upper[3] = {90, 108, 127.5};
  double emit[3] = {868008.225421772, 604533.07196845, 27253844.9047775};
  
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,iteration);
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_%i",sp->GetName(),(int)energy[en]));
  csp->SetTitle(Form("Sum bgsub xtal%i %.1fkeV",xtal,energy[en]));
  csp->Add(bg,-1);
  TH1D *ncsp = (TH1D *)NegRemove(csp);
  
  csp->GetXaxis()->SetRangeUser(lower[en]-20,upper[en]+20);
  new TCanvas;
  csp->Draw("hist");
  ncsp->SetLineColor(kRed);
  ncsp->Draw("hist same");

  TCanvas *c1 = new TCanvas;
  c1->SetTitle(Form("xtal%i bgsub10",xtal));
  TCanvas *c2 = new TCanvas;
  c2->SetTitle(Form("xtal%i bgsub10 negremove",xtal));
  c1->Divide(2,3);
  c2->Divide(2,3);
  std::vector<TF1 *> gpeak_bg;
  std::vector<TF1 *>gpeak_neg;
  for(int i=0;i<6;i++){
    double xlower = lower[en] - 0.5*i;
    double xupper = upper[en] + 0.5*i;
    csp = (TH1D *)sp->Clone(Form("%s_bgsub_%i",sp->GetName(),i));
    csp->Add(bg,-1);
    csp->SetTitle(Form("xtal%i bgsub10 [%.1f, %.1f]",xtal, xlower, xupper));
    csp->GetListOfFunctions()->Clear();
    csp->GetXaxis()->SetRangeUser(lower[en]-20, upper[en]+20);
    c1->cd(i+1);
    csp->Draw();
    if(en<2){
      gpeak_bg.push_back(GausFit(csp,xlower,xupper));
    }else{
      gpeak_bg.push_back(GausFit(csp,xlower,xupper));
    }
    c1->Update();    

    TH1D *ncsp = (TH1D *)NegRemove(csp);
    ncsp->SetTitle(Form("xtal%i bgsub10 negremove [%.1f, %.1f]",xtal, xlower, xupper));
    ncsp->GetListOfFunctions()->Clear();
    ncsp->GetXaxis()->SetRangeUser(lower[en]-20, upper[en]+20);
    c2->cd(i+1);
    ncsp->Draw();
    if(en<2){
      gpeak_neg.push_back(GausFit(ncsp,xlower,xupper));
    }else{
      gpeak_neg.push_back(GausFit(ncsp,xlower,xupper));
    }
    c2->Update();
  }

  std::cout<<"lower"<<"\t"<<"upper"<<"\t"
           <<"bgsub"<<"\t"<<"eff_bgsub""\t"
           <<"NegRemo"<<"\t"<<"eff_negremo""\t"
           <<"eff_dif(%)"<<std::endl;
  for(int i=0;i<gpeak_bg.size();i++){
    double xlower = lower[en] - 0.5*i;
    double xupper = upper[en] + 0.5*i;
    double sum_bg, sum_neg;
    if(en<2){
      sum_bg  = ((GGaus *)gpeak_bg[i])->GetSum();
      sum_neg = ((GGaus *)gpeak_neg[i])->GetSum();
    }else{
      sum_bg  = ((GPeak *)gpeak_bg[i])->GetSum();
      sum_neg = ((GPeak *)gpeak_neg[i])->GetSum();
    }
    double eff_bg = sum_bg/emit[en];
    double eff_neg = sum_neg/emit[en];
    double dif = fabs(eff_bg-eff_neg)/eff_bg*100;
    std::cout << xlower  << "\t"
              << xupper  << "\t"
              << sum_bg  << "\t"
              << eff_bg  << "\t"
              << sum_neg << "\t"
              << eff_neg << "\t"
              << dif     <<std::endl;
  }

  return;

}

//=================================================== SRMPeakFit ========================================//
//Fit all avaliable peaks from SRM
//Return sum from fitting as double *
//Return can be input for "SRMEffFit(double *det)"
double *SRMPeakFit(int xtal=0, int iteration=10, double *lower=0, double *upper=0, Option_t *opt = "bg"){
  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,iteration);
  
  TCanvas *c = new TCanvas;
  c->Divide(5,2);
  c->cd(1);
  sp->Draw();
  bg->Draw("same"); 

  int size = 9;
  double energy[9] = {105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower_raw[9]  = {102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double upper_raw[9] = {109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};
  
  if(lower==nullptr || upper==nullptr){
    lower = lower_raw;
    upper = upper_raw;
  }
  
  std::vector<TF1*> gpeak;
  std:vector<double> intg;
  for(int i=0;i<size;i++){
    TH1D *csp = (TH1D *)sp->Clone(Form("%s_%i",sp->GetName(),(int)energy[i]));
    csp->SetTitle(Form("Sum bgsub = %i xtal%i %.1fkeV",iteration,xtal,energy[i]));
    csp->Add(bg,-1);
    csp->GetListOfFunctions()->Clear();
    if(sopt.Contains("bg")){
      csp->GetXaxis()->SetRangeUser(lower[i]-20, upper[i]+20);
      intg.push_back(csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
      csp->GetXaxis()->FindBin(upper[i])));
      c->cd(i+2);
      csp->Draw();
      if(i<1){
        gpeak.push_back(GausFit(csp,lower[i],upper[i]));
      }else{
        gpeak.push_back(PhotoPeakFit(csp,lower[i],upper[i]));
      }
    }
    if(sopt.Contains("neg")){
      TH1D *ncsp = (TH1D *)NegRemove(csp);
      ncsp->GetXaxis()->SetRangeUser(lower[i]-20, upper[i]+20);
      intg.push_back(csp->Integral(ncsp->GetXaxis()->FindBin(lower[i]),
                     ncsp->GetXaxis()->FindBin(upper[i])));
      c->cd(i+2);
      ncsp->Draw();
      if(i<2){
        gpeak.push_back(GausFit(ncsp,lower[i],upper[i]));
      }else{
        gpeak.push_back(PhotoPeakFit(ncsp,lower[i],upper[i]));
      }
    }
    c->Update();
  }

  double *output = energy;
  std::cout<<"iteration = "<<iteration<<"\t"<<"opt = "<<sopt.Data()<<std::endl;
  std::cout<< "energy" << "\t";
  for(int i=0;i<size;i++){
    std::cout << energy[i] << "\t";
  }
  std::cout<<std::endl;
  std::cout<< "integ" << "\t";
  for(int i=0;i<size;i++){
    std::cout << intg[i] << "\t";
  }
  std::cout<<std::endl;
  std::cout<< "sum" << "\t";
  for(int i=0;i<size;i++){
    if(i<1){
      output[i] = ((GGaus *)gpeak[i])->GetSum();
      std::cout << ((GGaus *)gpeak[i])->GetSum() << "\t";
    }else{
      output[i] = ((GPeak *)gpeak[i])->GetSum();
      std::cout << ((GPeak *)gpeak[i])->GetSum() << "\t";
    }
  }
  std::cout<<std::endl; 

  return output;
}


//new Input: to change a specific peak fitting range to change the sum;
//Input en = the number of energy:
//0 = 86keV;
//1 = 105keV;
//2 = 123keV;
double *SRMPeakFit_Range(int en=0, double *fitrange=0, int xtal = 0, int iteration=10, Option_t *opt = "bg"){
  
  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");

  double lower[9]  = {102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double upper[9] = {109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};
  
  lower[en] = fitrange[0];
  upper[en] = fitrange[1];

  return SRMPeakFit(xtal, iteration, lower, upper, opt);
}



//============================================== SRMEffFit ============================================//
//Draw SRM efficiency plot
//Retunr efficiency function TF1 *
TGraph *SRMEffFit(double *det=0, Option_t *opt = "bg"){

  int size = 9;
  double energy[size];
  double energy_err[size];
  double eff[size];
  double eff_err[size];

  double det_bg[9] = {1646.48, 64740.7, 7067.47, 2591.63, 8766.86, 4528.81, 3466.57, 5932.3,9714.6};
  double det_neg[9] = {1299.4, 64731.8, 7055.4, 2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double energy_raw[9] = {105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double emit[9] = {604533.072, 27253844.9, 4618975.333, 3305907.141, 13414247.94, 8144646.136, 6970997.902, 12083063.03, 23299674.34}; 
  double emit_err[9] = {2.120788287, 1.313796825, 1.202523222, 1.202523222, 1.202523222, 1.255413119, 1.376975707, 1.255413119, 1.155881524};
  double sys_err = 1;

  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  if(det==nullptr){
    if(sopt.Contains("bg")) {
      det = det_bg;
    }
    if(sopt.Contains("neg")) {
      det = det_neg;
    }
  }
  
  for(int i=0;i<size;i++){
    energy[i] = energy_raw[i];
    energy_err[i] = 0;
    eff[i] = det[i]/emit[i];
    eff_err[i] = sqrt(1./det[i] + pow((sys_err/100),2)); //detecting error% + system error
    eff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(eff_err[i],2)); // detecting error% + emitting error%
    eff_err[i] = eff[i] * eff_err[i]; // absolute eff error
  }

  TGraphErrors *gr = new TGraphErrors(size,energy,eff,energy_err,eff_err);
  gr->SetTitle(Form("SRM Sum %s",sopt.Data()));
  TF1 *fx = new TF1("fxSRM", GammaEff, 0,4000,4);
  gr->Fit(fx);
  new TCanvas;
  gr->Draw("A*");

  std::cout << "Energy"    << "\t"
    << "Observed"  << "\t"
    << "Calculated"<< "\t"
    << "dif%"      << std::endl;
  for(int i=0;i<size;i++){
    double cal = fx->Eval(energy[i]);
    double dif = fabs(eff[i] - cal)/cal;
    std::cout << energy[i]<< "\t"
      << eff[i]   << "\t"
      << cal      << "\t"
      << dif      << std::endl;
  }

  return gr; 
}


//======================================= SRMEffFit_Range ============================//
//change one peak fit range to check how the efficiency changes;
//Input en = peak number
//0 = 105keV;
//1 = 123keV
void SRMEffFit_Range(int en=0, int xtal=0, int iteration=10){
  double lower[2] = {103.5, 118.5};
  double upper[2] = {108, 127.5};
 
  int size = 9;
  double energy[size];
  double energy_err[size];
  double bgeff[size];
  double bgeff_err[size];
  double negeff[size];
  double negeff_err[size];
 
  double det_bg[9] = {1646.48, 64740.7, 7067.47, 2591.63, 8766.86, 4528.81, 3466.57, 5932.3,9714.6};
  double det_neg[9] = {1299.4, 64731.8, 7055.4, 2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double energy_raw[9] = {105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double emit[9] = {604533.072, 27253844.9, 4618975.333, 3305907.141, 13414247.94, 8144646.136, 6970997.902, 12083063.03, 23299674.34};
  double emit_err[9] = {2.120788287, 1.313796825, 1.202523222, 1.202523222, 1.202523222, 1.255413119, 1.376975707, 1.255413119, 1.155881524};
  double sys_err = 1; 
  
  for(int i=0;i<size;i++){
    energy[i] = energy_raw[i];
    energy_err[i] = 0;
    bgeff[i] = det_bg[i]/emit[i];
    bgeff_err[i] = sqrt(1./det_bg[i] + pow((sys_err/100),2)); //detecting error% + system error
    bgeff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(bgeff_err[i],2)); // detecting error% + emitting error%
    bgeff_err[i] = bgeff[i] * bgeff_err[i]; // absolute eff error
    negeff[i] = det_neg[i]/emit[i];
    negeff_err[i] = sqrt(1./det_neg[i] + pow((sys_err/100),2)); //detecting error% + system error
    negeff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(negeff_err[i],2)); // detecting error% + emitting error%
    negeff_err[i] = negeff[i] * negeff_err[i]; // absolute eff error
  }
  



  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,iteration);
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_%i",sp->GetName(),(int)energy[en]));
  csp->SetTitle(Form("Sum bgsub xtal%i %.1fkeV",xtal,energy[en]));
  csp->Add(bg,-1);
  TH1D *ncsp = (TH1D *)NegRemove(csp);

  std::vector<double> sum_bg;
  std::vector<double> sum_neg;
  for(int i=0;i<6;i++){
    double xlower = lower[en] - 0.5*i;
    double xupper = upper[en] + 0.5*i;
    csp = (TH1D *)sp->Clone(Form("%s_bgsub_%i",sp->GetName(),i));
    csp->Add(bg,-1);
    csp->SetTitle(Form("xtal%i bgsub10 [%.1f, %.1f]",xtal, xlower, xupper));
    csp->GetListOfFunctions()->Clear();
    csp->GetXaxis()->SetRangeUser(lower[en]-20, upper[en]+20);
    if(en<2){
      sum_bg.push_back(GausFit(csp,xlower,xupper)->GetSum());
    }else{
      sum_bg.push_back(GausFit(csp,xlower,xupper)->GetSum());
    }

    TH1D *ncsp = (TH1D *)NegRemove(csp);
    ncsp->SetTitle(Form("xtal%i bgsub10 negremove [%.1f, %.1f]",xtal, xlower, xupper));
    ncsp->GetListOfFunctions()->Clear();
    ncsp->GetXaxis()->SetRangeUser(lower[en]-20, upper[en]+20);
    if(en<2){
      sum_neg.push_back(GausFit(ncsp,xlower,xupper)->GetSum());
    }else{
      sum_neg.push_back(GausFit(ncsp,xlower,xupper)->GetSum());
    }
  }

  TCanvas *c = new TCanvas;
  c->Divide(2,6);
  TF1 *fx[6];
  TF1 *fx1[6];
  for(int i=0;i<sum_bg.size();i++){
    double xlower = lower[en] - 0.5*i;
    double xupper = upper[en] + 0.5*i;
    bgeff[en] = sum_bg[i]/emit[en];
    bgeff_err[en] = 0;
    bgeff_err[en] = sqrt(1./sum_bg[i] + pow((sys_err/100),2)); 
    bgeff_err[en] = sqrt(pow((emit_err[en]/100),2) + pow(bgeff_err[en],2)); 
    bgeff_err[en] = bgeff[en] * bgeff_err[en]; 
    negeff[en] = sum_neg[i]/emit[en];
    negeff_err[en] = 0;   
    negeff_err[en] = sqrt(1./sum_neg[i] + pow((sys_err/100),2)); 
    negeff_err[en] = sqrt(pow((emit_err[en]/100),2) + pow(negeff_err[en],2)); 
    negeff_err[en] = negeff[en] * negeff_err[en];

    TGraphErrors *gr = new TGraphErrors(size,energy,bgeff,energy_err,bgeff_err);
    gr->SetName(Form("SRM_xtal%i_bg_%i",xtal,i));
    gr->SetTitle(Form("SRM xtal%i bgsub %.1f [%.1f,%.1f]",xtal,energy[en],xlower,xupper));
    fx[i] = new TF1(Form("fxSRM_bg_%i",i), GammaEff, 0,4000,4);
    gr->Fit(fx[i]); 
    TGraphErrors *gr1 = new TGraphErrors(size,energy,negeff,energy_err,negeff_err);
    gr1->SetName(Form("SRM_xtal%i_neg_%i",xtal,i));
    gr1->SetTitle(Form("SRM xtal%i negremove %.1f [%.1f,%.1f]",xtal,energy[en],xlower,xupper));
    fx1[i] = new TF1(Form("fxSRM_neg_%i",i), GammaEff, 0,4000,4);
    gr1->Fit(fx1[i]); 
    c->cd(2*i+1);
    gr->Draw("A*");
    c->cd(2*i+2);
    gr1->Draw("A*");
    c->Update();
    

  }
  

  std::cout << "lower" << "\t"  << "upper"  << "\t"
            << "sum_bg" << "\t" << "observed_bg" << "\t" << "calculate_bg" << "\t" << "dif_bg" << "\t"
            << "sum_neg"<< "\t" << "observed_neg"<< "\t" << "calculate_neg"<< "\t" << "dif_neg"<<std::endl;
  for(int i=0;i<6;i++){
    double xlower = lower[en] - 0.5*i;
    double xupper = upper[en] + 0.5*i;
    double cal = fx[i]->Eval(energy[en]);
    double dif = fabs(bgeff[en] - cal)/cal*100; 
    double cal1 = fx1[i]->Eval(energy[en]);
    double dif1 = fabs(negeff[en] - cal1)/cal1*100; 
    std::cout << xlower << "\t" << xupper << "\t"
              << sum_bg[i] << "\t" << bgeff[en] << "\t" << cal << "\t" << dif << "\t"
              << sum_neg[i] << "\t" << negeff[en] << "\t" << cal1 << "\t" << dif1 << std::endl;
  }
  

  return; 
}



