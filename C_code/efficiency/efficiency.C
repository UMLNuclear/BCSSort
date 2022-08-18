
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


//=========================== SRMFit =========================//
TF1 *SRMFit(double *det=0,Option_t *opt="xtal0",double sys_err=3){
  
  int size = 9;
  double energy_raw[9] = {105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double emit[9] = {604533.072,27253844.9,4618975.333,3305907.141,13414247.94,8144646.136,6970997.902,12083063.03,23299674.34};
  double emit_err[9]={2.120788287,1.313796825,1.202523222,1.202523222,1.202523222,1.255413119,1.376975707,1.255413119,1.155881524};
  
  double det_xtal0[9] = {1299.4,  64731.8, 7055.4,  2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double det_xtal1[9] = {1605.86, 66260.7, 6950.57, 2507.2,  8773.72, 4388.24, 3392.1,  5809.28, 9533.45};
  double det_xtal2[9] = {1597.92, 68887.2, 7595.84, 2567.54, 9331.19, 4616.28, 3466.58, 6182.31, 10211.9};
  double det_xtal3[9] = {1622.82, 70125.3, 7420.91, 2534.19, 9130.54, 4757.94, 3576.49, 6080.55, 10200.7};
  double det_sum0[9]  = {5819.355924, 269779.004785,  29019.8, 10132,   35743.3, 18411.3, 13940.1, 24092.4, 39726.3};
  double det_addb0[9] = {4782,    266227,  31057.5, 12881.2, 49444.5, 25953.6, 19156.5, 33564.6, 58687.9};
  

  double raw_xtal0[9] = {3192.87, 66513.6, 7869.48, 2769.16, 8713.03, 4735.37, 3625.85, 6082.4 , 9894.35};
  double raw_xtal1[9] = {3018.44, 67977.6, 7835.49, 2684.95, 8900.74, 4532.02, 3549.24, 6089.47, 9660.47};
  double raw_xtal2[9] = {3011.62, 71374.1, 8455.41, 2760.45, 9129.08, 4768.5 , 3732.44, 6420.5 , 10330.2};
  double raw_xtal3[9] = {3472.39, 71825.7, 8588.53, 2824.01, 8993.74, 4967.8 , 3840.89, 6336.46, 10333.1};
  double raw_sum0[9]  = {12595.1, 277340 , 32470.2, 10985.6, 35670.3, 18962.7, 14681.1, 24887.7, 40121};
  double raw_addb0[9] = {11770.4, 272361 , 33883.9, 14020.1, 47985.8, 26672.8, 19628.2, 35769.8, 59905.2};

  double energy[size];
  double energy_err[size];
  double eff[size];
  double eff_err[size];


  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  if(det==nullptr){
    if(sopt.Contains("x") && sopt.Contains("0")){
      det = det_xtal0;
    }
    if(sopt.Contains("x") && sopt.Contains("1")){
      det = det_xtal1;
    }
    if(sopt.Contains("x") && sopt.Contains("2")){
      det = det_xtal2;
    }
    if(sopt.Contains("x") && sopt.Contains("3")){
      det = det_xtal0;
    }
    if(sopt.Contains("sum")){
      det = det_sum0;
    }
    if(sopt.Contains("add")){
      det = det_addb0;
    }
    if(sopt.Contains("raw")){
      if(sopt.Contains("x") && sopt.Contains("0")){
        det = raw_xtal0;
      }
      if(sopt.Contains("x") && sopt.Contains("1")){
        det = raw_xtal1;
      }
      if(sopt.Contains("x") && sopt.Contains("2")){
        det = raw_xtal2;
      }
      if(sopt.Contains("x") && sopt.Contains("3")){
        det = raw_xtal0;
      }
      if(sopt.Contains("sum")){
        det = raw_sum0;
      }
      if(sopt.Contains("add")){
        det = raw_addb0;
      }
    }
  }

  for(int i=0;i<size;i++){
    //==============================//
    if(energy[i]>500) sys_err = 1;
    //==============================//
    energy[i] = energy_raw[i];
    energy_err[i] = 0;
    eff[i] = det[i]/emit[i];
    eff_err[i] = sqrt(1./det[i] + pow((sys_err/100),2)); //detecting error% + system error
    eff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(eff_err[i],2)); // detecting error% + emitting error%
    eff_err[i] = eff[i] * eff_err[i]; // absolute eff error
  }
  TGraphErrors *gr = new TGraphErrors(size,energy,eff,energy_err,eff_err);
  gr->SetTitle(Form("SRM %s",sopt.Data()));
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
    double dif = fabs(eff[i] - cal)/cal*100;
    std::cout << energy[i]<< "\t"
      << eff[i]   << "\t"
      << cal      << "\t"
      << dif      << std::endl;
  }

  return fx;
}

//================================== Efficicency Fit with SRM and 56Co =========================//
TF1 *EffFit(double *detsrm=0, double *detco=0, Option_t* opt="addb", double sys_err=3){
  
  int ssrm = 9;
  int sco  = 12;
  double energy_srm[9] = {105.3,   123.1,   247.7,    591.8,    723.3,    873.2,    996.3,    1004.7,   1274.5};
  double energy_co[12]  =  {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  double inten_co[12]   =  {100,     14.13,   2.238,    66.07,    4.255,    15.48,    3.029,    7.77,     16.96,    3.13,     7.62,     1.78};
  double emit_srm[9] = {604533.072,27253844.9,4618975.333,3305907.141,13414247.94,8144646.136,6970997.902,12083063.03,23299674.34};
  double emit_srm_err[9]={2.120788287,1.313796825,1.202523222,1.202523222,1.202523222,1.255413119,1.376975707,1.255413119,1.155881524};
  //double sys_err = 1;
  //double sys_err = 5;
  
  double detco_xtal0[12] = {598134,  72222.5, 10085.2, 299019,  17875.4, 50083.9, 8896.92, 23142.8, 40021.8, 5904.64, 14343.4, 3308.78};
  double detco_xtal1[12] = {598426,  71888.4, 10343.1, 299037,  18011,   49794,   8832.74, 22602.7, 40070.8, 5836.12, 14221.1, 3408.68};
  double detco_xtal2[12] = {643877,  77630.2, 11079.1, 322008,  18813.3, 54080.9, 9606.03, 24742.2, 42984.1, 6410.33, 15709.9, 3680.42};
  double detco_xtal3[12] = {637837,  76626.2, 11206,   315781,  18947.2, 53740.3, 9703.44, 24323.7, 41847.8, 6210.19, 14835.1, 3570.83};
  double detco_sum0[12]  = {2478500, 298517,  42793.8, 1237340, 73846.3, 207777,  37093.4, 94844.1, 164886 , 24394.1, 59050.7, 13939.7};
  double detco_addb0[12] = {3430810, 421095,  61831.1, 1768570, 111909,  308067,  56042.4, 143321 , 258248 , 39177.2, 95393.5, 22781.6};

  double rawco_xtal0[12] = {600808 , 73923.2, 10927.1, 299940 , 18492.7, 50977.1, 9258.85, 23400.1, 40778.5, 6061.62, 14691.3, 3368.2};
  double rawco_xtal1[12] = {600588 , 73165.2, 11127.3, 299232 , 18530.5, 50657.5, 9126.58, 22820.3, 40582.8, 6024.99, 14448.2, 3480.06};
  double rawco_xtal2[12] = {650757 , 79321.6, 12008.3, 324115 , 19623.4, 55592.9, 9882.23, 25253.7, 44050.8, 6583.87, 16202.3, 3773.65};
  double rawco_xtal3[12] = {641503 , 78234.8, 12070.8, 318744 , 19899.2, 56366.3, 10206.5, 25240.9, 43825.7, 6447.5 , 15260.5, 3786.13};
  double rawco_sum0[12]  = {2476010, 304041,  46019.4, 1241830, 76473.3, 216189 , 38619.3, 96651.4, 167466,  25126  , 60558.8, 14455.8};
  double rawco_addb0[12] = {3474940, 426750,  64703.1, 1778550, 116105 , 326390 , 58778.5, 148734 , 272684,  42583.1, 103551 , 24698.4};
  
  double detsrm_xtal0[9] = {1299.4,  64731.8, 7055.4,  2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double detsrm_xtal1[9] = {1605.86, 66260.7, 6950.57, 2507.2,  8773.72, 4388.24, 3392.1,  5809.28, 9533.45};
  double detsrm_xtal2[9] = {1597.92, 68887.2, 7595.84, 2567.54, 9331.19, 4616.28, 3466.58, 6182.31, 10211.9};
  double detsrm_xtal3[9] = {1622.82, 70125.3, 7420.91, 2534.19, 9130.54, 4757.94, 3576.49, 6080.55, 10200.7};
  //double detsrm_sum0[9]  = {5819.355924, 269779.004785, 29019.8, 10132,   35743.3, 18411.3, 13940.1, 24092.4, 39726.3};
  double detsrm_sum0[9]  = {5819.355924, 269779.004785, 29095.195959, 10132,   35743.3, 18411.3, 13940.1, 24092.4, 39726.3};
  double detsrm_addb0[9] = {5536.5, 265875.91, 31057.5, 12881.2, 49444.5, 25953.6, 19156.5, 33564.6, 58687.9};
  //double detsrm_addb0[9] = {4782,    266227,  31057.5, 12881.2, 49444.5, 25953.6, 19156.5, 33564.6, 58687.9};

  double rawsrm_xtal0[9] = {3192.87, 66513.6, 7869.48, 2769.16, 8713.03, 4735.37, 3625.85, 6082.4 , 9894.35};
  double rawsrm_xtal1[9] = {3018.44, 67977.6, 7835.49, 2684.95, 8900.74, 4532.02, 3549.24, 6089.47, 9660.47};
  double rawsrm_xtal2[9] = {3011.62, 71374.1, 8455.41, 2760.45, 9129.08, 4768.5 , 3732.44, 6420.5 , 10330.2};
  double rawsrm_xtal3[9] = {3472.39, 71825.7, 8588.53, 2824.01, 8993.74, 4967.8 , 3840.89, 6336.46, 10333.1};
  double rawsrm_sum0[9]  = {12595.1, 277340 , 32470.2, 10985.6, 35670.3, 18962.7, 14681.1, 24887.7, 40121};
  double rawsrm_addb0[9] = {11770.4, 272361 , 33883.9, 14020.1, 47985.8, 26672.8, 19628.2, 35769.8, 59905.2};
  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  if(detco==nullptr){
    if(sopt.Contains("x") && sopt.Contains("0")){
      detco  = detco_xtal0;
      detsrm = detsrm_xtal0;
    }
    if(sopt.Contains("x") && sopt.Contains("1")){
      detco  = detco_xtal1;
      detsrm = detsrm_xtal1;
    }
    if(sopt.Contains("x") && sopt.Contains("2")){
      detco  = detco_xtal2;
      detsrm = detsrm_xtal2;
    }
    if(sopt.Contains("x") && sopt.Contains("3")){
      detco  = detco_xtal3;
      detsrm = detsrm_xtal3;
    }
    if(sopt.Contains("sum")){
      detco  = detco_sum0;
      detsrm = detsrm_sum0;
    }
    if(sopt.Contains("add")){
      detco  = detco_addb0;
      detsrm = detsrm_addb0;
    }
    if(sopt.Contains("raw")){
      if(sopt.Contains("x") && sopt.Contains("0")){
        detco  = rawco_xtal0;
        detsrm = rawsrm_xtal0;
      }
      if(sopt.Contains("x") && sopt.Contains("1")){
        detco  = rawco_xtal1;
        detsrm = rawsrm_xtal1;
      }
      if(sopt.Contains("x") && sopt.Contains("2")){
        detco  = rawco_xtal2;
        detsrm = rawsrm_xtal2;
      }
      if(sopt.Contains("x") && sopt.Contains("3")){
        detco  = rawco_xtal3;
        detsrm = rawsrm_xtal3;
      }
      if(sopt.Contains("sum")){
        detco  = rawco_sum0;
        detsrm = rawsrm_sum0;
      }
      if(sopt.Contains("add")){
        detco  = rawco_addb0;
        detsrm = rawsrm_addb0;
      }
    }
  }

  TF1 *feff = SRMFit(detsrm,opt,sys_err);
  double emit_co = 0;
  for(int i=0;i<4;i++){
    double eff_temp = feff->Eval(energy_co[i]);
    emit_co += detco[i]/eff_temp/inten_co[i];
  }
  emit_co = emit_co/4.;    
 
  int size = ssrm + sco;
  double emit[size]; 
  double emit_err[size]; 
  double energy[size]; 
  double energy_err[size]; 
  double eff[size]; 
  double eff_err[size]; 
  double det[size];

  for(int i=0;i<size;i++){
    if(i<ssrm){
      energy[i] = energy_srm[i];
      det[i] = detsrm[i];
      emit[i] = emit_srm[i];
      emit_err[i] = emit_srm_err[i];
    }else{
      energy[i] = energy_co[i-ssrm];
      det[i] = detco[i-ssrm];
      emit[i] = emit_co*inten_co[i-ssrm];
      emit_err[i] = 1.5; //estimate the emit error for 56Co is about 1.5%.
    }
  }
  
  for(int i=0;i<size;i++){
    //===========//
    if(energy[i]>500) sys_err = 1;
    //==========//
    energy_err[i] = 0;
    eff[i] = det[i]/emit[i];
    eff_err[i] = sqrt(1./det[i] + pow((sys_err/100),2)); //detecting error% + system error
    eff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(eff_err[i],2)); // detecting error% + emitting error%
    eff_err[i] = eff[i] * eff_err[i]; // absolute eff error
  }

  TGraphErrors *gr = new TGraphErrors(size,energy,eff,energy_err,eff_err);
  gr->SetTitle(Form("SRM+56Co %s",sopt.Data()));
  TF1 *fx = new TF1("fx", GammaEff, 0,4000,4);
  gr->Fit(fx);
  new TCanvas;
  gr->Draw("A*");

  std::cout << "Energy"    << "\t"
    << "Observed"  << "\t"
    << "Calculated"<< "\t"
    << "dif%"      << std::endl;
  for(int i=0;i<size;i++){
    double cal = fx->Eval(energy[i]);
    double dif = fabs(eff[i] - cal)/cal*100;
    std::cout << energy[i]<< "\t"
      << eff[i]   << "\t"
      << cal      << "\t"
      << dif      << std::endl;
  }
  
  return fx;

}



//=============================================================================================================//



TF1 *SRMFit_no105(double *det=0, Option_t *opt="xtal0"){
  
  int size = 8;
  double energy_raw[8] = {123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double emit[8] = {27253844.9,4618975.333,3305907.141,13414247.94,8144646.136,6970997.902,12083063.03,23299674.34};
  double emit_err[8]={1.313796825,1.202523222,1.202523222,1.202523222,1.255413119,1.376975707,1.255413119,1.155881524};
  double sys_err = 1;
  
  double det_xtal0[8] = {64731.8, 7055.4,  2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double det_xtal1[8] = {66260.7, 6950.57, 2507.2,  8773.72, 4388.24, 3392.1,  5809.28, 9533.45};
  double det_xtal2[8] = {68887.2, 7595.84, 2567.54, 9331.19, 4616.28, 3466.58, 6182.31, 10211.9};
  double det_xtal3[8] = {70125.3, 7420.91, 2534.19, 9130.54, 4757.94, 3576.49, 6080.55, 10200.7};
  double det_sum0[8]  = {270045,  29019.8, 10132,   35743.3, 18411.3, 13940.1, 24092.4, 39726.3};
  double det_addb0[8] = {266227,  31057.5, 12881.2, 49444.5, 25953.6, 19156.5, 33564.6, 58687.9};

  double raw_xtal0[8] = {66513.6, 7869.48, 2769.16, 8713.03, 4735.37, 3625.85, 6082.4 , 9894.35};
  double raw_xtal1[8] = {67977.6, 7835.49, 2684.95, 8900.74, 4532.02, 3549.24, 6089.47, 9660.47};
  double raw_xtal2[8] = {71374.1, 8455.41, 2760.45, 9129.08, 4768.5 , 3732.44, 6420.5 , 10330.2};
  double raw_xtal3[8] = {71825.7, 8588.53, 2824.01, 8993.74, 4967.8 , 3840.89, 6336.46, 10333.1};
  double raw_sum0[8]  = {277340 , 32470.2, 10985.6, 35670.3, 18962.7, 14681.1, 24887.7, 40121};
  double raw_addb0[8] = {272361 , 33883.9, 14020.1, 47985.8, 26672.8, 19628.2, 35769.8, 59905.2};

  double energy[size];
  double energy_err[size];
  double eff[size];
  double eff_err[size];


  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  if(det==nullptr){
    if(sopt.Contains("x") && sopt.Contains("0")){
      det = det_xtal0;
    }
    if(sopt.Contains("x") && sopt.Contains("1")){
      det = det_xtal1;
    }
    if(sopt.Contains("x") && sopt.Contains("2")){
      det = det_xtal2;
    }
    if(sopt.Contains("x") && sopt.Contains("3")){
      det = det_xtal0;
    }
    if(sopt.Contains("sum")){
      det = det_sum0;
    }
    if(sopt.Contains("add")){
      det = det_addb0;
    }
    if(sopt.Contains("raw")){
      if(sopt.Contains("x") && sopt.Contains("0")){
        det = raw_xtal0;
      }
      if(sopt.Contains("x") && sopt.Contains("1")){
        det = raw_xtal1;
      }
      if(sopt.Contains("x") && sopt.Contains("2")){
        det = raw_xtal2;
      }
      if(sopt.Contains("x") && sopt.Contains("3")){
        det = raw_xtal0;
      }
      if(sopt.Contains("sum")){
        det = raw_sum0;
      }
      if(sopt.Contains("add")){
        det = raw_addb0;
      }
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
  gr->SetTitle(Form("SRM %s",sopt.Data()));
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
    double dif = fabs(eff[i] - cal)/cal*100;
    std::cout << energy[i]<< "\t"
      << eff[i]   << "\t"
      << cal      << "\t"
      << dif      << std::endl;
  }

  return fx;
}











TF1 *EffFit_no105(double *detco=0, double *detsrm=0, Option_t* opt="xtal0"){

  
  int ssrm = 8;
  int sco  = 12;
  double energy_srm[8] = {123.1,   247.7,    591.8,    723.3,    873.2,    996.3,    1004.7,   1274.5};
  double energy_co[12]  =  {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  double inten_co[12]   =  {100,     14.13,   2.238,    66.07,    4.255,    15.48,    3.029,    7.77,     16.96,    3.13,     7.62,     1.78};
  double emit_srm[8] = {27253844.9,4618975.333,3305907.141,13414247.94,8144646.136,6970997.902,12083063.03,23299674.34};
  double emit_srm_err[8]={1.313796825,1.202523222,1.202523222,1.202523222,1.255413119,1.376975707,1.255413119,1.155881524};
  double sys_err = 1;
  
  double detco_xtal0[12] = {598134,  72222.5, 10085.2, 299019,  17875.4, 50083.9, 8896.92, 23142.8, 40021.8, 5904.64, 14343.4, 3308.78};
  double detco_xtal1[12] = {598426,  71888.4, 10343.1, 299037,  18011,   49794,   8832.74, 22602.7, 40070.8, 5836.12, 14221.1, 3408.68};
  double detco_xtal2[12] = {643877,  77630.2, 11079.1, 322008,  18813.3, 54080.9, 9606.03, 24742.2, 42984.1, 6410.33, 15709.9, 3680.42};
  double detco_xtal3[12] = {637837,  76626.2, 11206,   315781,  18947.2, 53740.3, 9703.44, 24323.7, 41847.8, 6210.19, 14835.1, 3570.83};
  double detco_sum0[12]  = {2478500, 298517,  42793.8, 1237340, 73846.3, 207777,  37093.4, 94844.1, 164886 , 24394.1, 59050.7, 13939.7};
  double detco_addb0[12] = {3430810, 421095,  61831.1, 1768570, 111909,  308067,  56042.4, 143321 , 258248 , 39177.2, 95393.5, 22781.6};

  double rawco_xtal0[12] = {600808 , 73923.2, 10927.1, 299940 , 18492.7, 50977.1, 9258.85, 23400.1, 40778.5, 6061.62, 14691.3, 3368.2};
  double rawco_xtal1[12] = {600588 , 73165.2, 11127.3, 299232 , 18530.5, 50657.5, 9126.58, 22820.3, 40582.8, 6024.99, 14448.2, 3480.06};
  double rawco_xtal2[12] = {650757 , 79321.6, 12008.3, 324115 , 19623.4, 55592.9, 9882.23, 25253.7, 44050.8, 6583.87, 16202.3, 3773.65};
  double rawco_xtal3[12] = {641503 , 78234.8, 12070.8, 318744 , 19899.2, 56366.3, 10206.5, 25240.9, 43825.7, 6447.5 , 15260.5, 3786.13};
  double rawco_sum0[12]  = {2476010, 304041,  46019.4, 1241830, 76473.3, 216189 , 38619.3, 96651.4, 167466,  25126  , 60558.8, 14455.8};
  double rawco_addb0[12] = {3474940, 426750,  64703.1, 1778550, 116105 , 326390 , 58778.5, 148734 , 272684,  42583.1, 103551 , 24698.4};
  
  double detsrm_xtal0[8] = {64731.8, 7055.4,  2539.84, 8550.37, 4580.83, 3450.42, 5937.54, 9770.72};
  double detsrm_xtal1[8] = {66260.7, 6950.57, 2507.2,  8773.72, 4388.24, 3392.1,  5809.28, 9533.45};
  double detsrm_xtal2[8] = {68887.2, 7595.84, 2567.54, 9331.19, 4616.28, 3466.58, 6182.31, 10211.9};
  double detsrm_xtal3[8] = {70125.3, 7420.91, 2534.19, 9130.54, 4757.94, 3576.49, 6080.55, 10200.7};
  double detsrm_sum0[8]  = {270045,  29019.8, 10132,   35743.3, 18411.3, 13940.1, 24092.4, 39726.3};
  double detsrm_addb0[8] = {266227,  31057.5, 12881.2, 49444.5, 25953.6, 19156.5, 33564.6, 58687.9};

  double rawsrm_xtal0[9] = {3192.87, 66513.6, 7869.48, 2769.16, 8713.03, 4735.37, 3625.85, 6082.4 , 9894.35};
  double rawsrm_xtal1[9] = {3018.44, 67977.6, 7835.49, 2684.95, 8900.74, 4532.02, 3549.24, 6089.47, 9660.47};
  double rawsrm_xtal2[9] = {3011.62, 71374.1, 8455.41, 2760.45, 9129.08, 4768.5 , 3732.44, 6420.5 , 10330.2};
  double rawsrm_xtal3[9] = {3472.39, 71825.7, 8588.53, 2824.01, 8993.74, 4967.8 , 3840.89, 6336.46, 10333.1};
  double rawsrm_sum0[9]  = {12595.1, 277340 , 32470.2, 10985.6, 35670.3, 18962.7, 14681.1, 24887.7, 40121};
  double rawsrm_addb0[9] = {11770.4, 272361 , 33883.9, 14020.1, 47985.8, 26672.8, 19628.2, 35769.8, 59905.2};
  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  if(detco==nullptr){
    if(sopt.Contains("x") && sopt.Contains("0")){
      detco  = detco_xtal0;
      detsrm = detsrm_xtal0;
    }
    if(sopt.Contains("x") && sopt.Contains("1")){
      detco  = detco_xtal1;
      detsrm = detsrm_xtal1;
    }
    if(sopt.Contains("x") && sopt.Contains("2")){
      detco  = detco_xtal2;
      detsrm = detsrm_xtal2;
    }
    if(sopt.Contains("x") && sopt.Contains("3")){
      detco  = detco_xtal3;
      detsrm = detsrm_xtal3;
    }
    if(sopt.Contains("sum")){
      detco  = detco_sum0;
      detsrm = detsrm_sum0;
    }
    if(sopt.Contains("add")){
      detco  = detco_addb0;
      detsrm = detsrm_addb0;
    }
    if(sopt.Contains("raw")){
      if(sopt.Contains("x") && sopt.Contains("0")){
        detco  = rawco_xtal0;
        detsrm = rawsrm_xtal0;
      }
      if(sopt.Contains("x") && sopt.Contains("1")){
        detco  = rawco_xtal1;
        detsrm = rawsrm_xtal1;
      }
      if(sopt.Contains("x") && sopt.Contains("2")){
        detco  = rawco_xtal2;
        detsrm = rawsrm_xtal2;
      }
      if(sopt.Contains("x") && sopt.Contains("3")){
        detco  = rawco_xtal3;
        detsrm = rawsrm_xtal3;
      }
      if(sopt.Contains("sum")){
        detco  = rawco_sum0;
        detsrm = rawsrm_sum0;
      }
      if(sopt.Contains("add")){
        detco  = rawco_addb0;
        detsrm = rawsrm_addb0;
      }
    }
  }
  
  TF1 *feff = SRMFit_no105(detsrm,opt);
  double emit_co = 0;
  for(int i=0;i<4;i++){
    double eff_temp = feff->Eval(energy_co[i]);
    emit_co += detco[i]/eff_temp/inten_co[i];
  }
  emit_co = emit_co/4.;    
 
  int size = ssrm + sco;
  double emit[size]; 
  double emit_err[size]; 
  double energy[size]; 
  double energy_err[size]; 
  double eff[size]; 
  double eff_err[size]; 
  double det[size];

  for(int i=0;i<size;i++){
    if(i<ssrm){
      energy[i] = energy_srm[i];
      det[i] = detsrm[i];
      emit[i] = emit_srm[i];
      emit_err[i] = emit_srm_err[i];
    }else{
      energy[i] = energy_co[i-ssrm];
      det[i] = detco[i-ssrm];
      emit[i] = emit_co*inten_co[i-ssrm];
      emit_err[i] = 1.5; //estimate the emit error for 56Co is about 1.5%.
    }
  }
  
  for(int i=0;i<size;i++){
    energy_err[i] = 0;
    eff[i] = det[i]/emit[i];
    eff_err[i] = sqrt(1./det[i] + pow((sys_err/100),2)); //detecting error% + system error
    eff_err[i] = sqrt(pow((emit_err[i]/100),2) + pow(eff_err[i],2)); // detecting error% + emitting error%
    eff_err[i] = eff[i] * eff_err[i]; // absolute eff error
  }

  TGraphErrors *gr = new TGraphErrors(size,energy,eff,energy_err,eff_err);
  gr->SetTitle(Form("SRM+56Co %s",sopt.Data()));
  TF1 *fx = new TF1("fx", GammaEff, 0,4000,4);
  gr->Fit(fx);
  new TCanvas;
  gr->Draw("A*");

  std::cout << "Energy"    << "\t"
    << "Observed"  << "\t"
    << "Calculated"<< "\t"
    << "dif%"      << std::endl;
  for(int i=0;i<size;i++){
    double cal = fx->Eval(energy[i]);
    double dif = fabs(eff[i] - cal)/cal*100;
    std::cout << energy[i]<< "\t"
      << eff[i]   << "\t"
      << cal      << "\t"
      << dif      << std::endl;
  }
  
  return fx;

}
