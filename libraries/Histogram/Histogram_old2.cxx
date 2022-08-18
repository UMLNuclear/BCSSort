
#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<map>
#include<string>

#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TCutG.h>

#include <Histogram.h>
#include <OutputManager.h>
#include <util.h>
#include <TChannel.h>
#include <TOFCorrection.h>
#include <DetHit.h>
#include <Implant.h>
#include <BCSint.h>
#include <ddaschannel.h>




//======================== Constructor & Destructor========================//


Histogram *Histogram::fHistogram = 0;

Histogram *Histogram::Get(){
  if(!fHistogram){
    fHistogram = new Histogram;
  }
  return fHistogram;
}

Histogram::Histogram(){};
Histogram::~Histogram(){};

//========================== Beta Sort ===============================//

void Histogram::BetaTOF(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();
  runnum = runnum.substr(0,4);
  std::vector<double> toflin_Na = TOFCorrection::Get()->ReadFile(std::stoi(runnum), "/home/zhu/packages/BCSSort/config/TOF/TOF_beta_Na_lin.txt");
  std::vector<double> tofoff_Na = TOFCorrection::Get()->ReadFile(std::stoi(runnum), "/home/zhu/packages/BCSSort/config/TOF/TOF_beta_Na_offset.txt");


  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Na = (TCutG *)cutf->Get("Na"); 
  TCutG *Ne = (TCutG *)cutf->Get("Ne"); 

  long x=0;
  long n = beta->GetEntries();

  for(x=0;x<n;x++){
    beta->GetEntry(x);
    FillHistogram("tof",500,11000,16000,fbeta->fImplant.fI2S);
    FillHistogram("tof_entry",n,0,n,x, 500,11000,16000,fbeta->fImplant.fI2S);
    if(Na->IsInside(fbeta->fImplant.fI2S, fbeta->fImplant.fPIN1E)){
      double tofNa = fbeta->fImplant.fI2S;
      FillHistogram("tof_Na_before", 500,11000,16000,tofNa);
      double tofNa1 = toflin_Na[0] + toflin_Na[1]*tofNa;
      FillHistogram("tof_Na_linear", 500,11000,16000,tofNa1);
      double tofNa0 = tofoff_Na[0] + tofNa;
      FillHistogram("tof_Na_offset", 500,11000,16000,tofNa0);
    }
    if(Ne->IsInside(fbeta->fImplant.fI2S, fbeta->fImplant.fPIN1E)){
      FillHistogram("tof_Ne", 900,11000,15500,fbeta->fImplant.fI2S);
    }
    if((x%20000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }

  }

  printf("   on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("beta_tof%s.root",runnum.c_str()));

}


void Histogram::Beta3DPID(){
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation1bestT/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  long x=0;
  long entries = beta->GetEntries();
  double tofbd = 8000;

  for(x=0;x<entries;x++){
    fbeta->Clear();
    beta->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    double tof = fimp.fI2S;
    double i2 = fimp.fI2S_I2N;
    double pin1e = fimp.fPIN1E;
    for(int d=0;d<fbeta->DecaySize();d++){
      double decaytime = fbeta->fDecay[d].fDecayTime;
      if(decaytime>500) continue;
      FillHistogram("PID_dt",800,tofbd,tofbd+8000,tof, 300,4500,7500,fimp.fPIN1E, 500,0,500,decaytime);
      for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
        double e = fbeta->fDecay[d].fGe[m].GetEnergy();
        if(decaytime>500) continue;
        if(e>=882 && e<=888){
          FillHistogram("PID_dt_885",800,tofbd,tofbd+8000,tof, 300,4500,7500,fimp.fPIN1E, 500,0,500,decaytime);
        }       
        if(e>=900 && e<=906){
          FillHistogram("PID_dt_BGR885",800,tofbd,tofbd+8000,tof, 300,4500,7500,fimp.fPIN1E, 500,0,500,decaytime);
        }       
        if(decaytime>300) continue;
        if(e>=146 && e<=152){
          FillHistogram("PID_dt_150",800,tofbd,tofbd+8000,tof, 300,4500,7500,fimp.fPIN1E, 500,0,500,decaytime);
        }       
        if(e>=156 && e<=162){
          FillHistogram("PID_dt_BGR150",800,tofbd,tofbd+8000,tof, 300,4500,7500,fimp.fPIN1E, 500,0,500,decaytime);
        }       
      }
    } 
    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,entries);
      fflush(stdout);
    }
  }
  printf("   on entry %lu / %lu   \n",x,entries);
  SaveHistograms(Form("beta_prompt_op%s.root",runnum.c_str()));

  return;
}



void Histogram::Beta150GateTOF(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation1bestT/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  long x=0;
  long entries = beta->GetEntries();
  int tofbd = 6000;
  double slope = -0.24190693;
  double offset = 9071.3105;
  double tof2;
  double partof2  = -0.257303; //tan(TOF = tan*i2)
  for(0;x<entries;x++){
    fbeta->Clear();
    beta->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    double tof = fimp.fI2S;
    double i2 = fimp.fI2S_I2N;
    double pin1e = fimp.fPIN1E;
    double jug = slope*tof + offset;
    for(int d=0;d<fbeta->DecaySize();d++){
      double decaytime = fbeta->fDecay[d].fDecayTime;
      if(decaytime>60) continue;
      for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
        double e = fbeta->fDecay[d].fGe[m].GetEnergy();
        if(e<10 || e>4000) continue;
        if(pin1e<jug){
          if(decaytime<30){
            FillHistogram("singles30ms_TOF_Ne",2000,tofbd,tofbd+20000,tof, 4000,0,4000,e);
          }
          if(e>=146 && e<=152){//150keV
            if(decaytime<30){
              FillHistogram("pid_30ms_150",   2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("dt_TOF_30ms_150",2000,tofbd,tofbd+20000,tof, 3000,0,30,decaytime);  
              if(i2>10000 && i2<16000){
                tof2 = tof - partof2*i2; 
                FillHistogram("pid_30ms_150_tofCor",   2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
                FillHistogram("dt_TOF_30ms_150_tofCor",2000,tofbd,tofbd+20000,tof2, 3000,0,30,decaytime);  
              }
            }
          }
          if(e>=154 && e<=160){//BGR_150keV
            if(decaytime<30){
              FillHistogram("pid_30ms_BGR150",   2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("dt_TOF_30ms_BGR150",2000,tofbd,tofbd+20000,tof, 3000,0,30,decaytime);  
              if(i2>10000 && i2<16000){
                tof2 = tof - partof2*i2; 
                FillHistogram("pid_30ms_BGR150_tofCor",   2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
                FillHistogram("dt_TOF_30ms_bgR150_tofCor",2000,tofbd,tofbd+20000,tof2, 3000,0,30,decaytime);  
              }
            }
          }
          if(e>=152 && e<=188){//BGR_150keV
            if(decaytime<30){
              FillHistogram("pid_30ms_extBGR150",   2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("dt_TOF_30ms_extBGR150",2000,tofbd,tofbd+20000,tof, 3000,0,30,decaytime);  
              if(i2>10000 && i2<16000){
                tof2 = tof - partof2*i2; 
                FillHistogram("pid_30ms_extBGR150_tofCor",   2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
                FillHistogram("dt_TOF_30ms_extBGR150_tofCor",2000,tofbd,tofbd+20000,tof2, 3000,0,30,decaytime);  
              }
            }
          }
        }else{
          if(decaytime<60){
            FillHistogram("singles60ms_TOF_Na",2000,tofbd,tofbd+20000,tof, 4000,0,4000,e);
          }
          if(e>=882 && e<=888){//885keV
            if(decaytime<60){
              FillHistogram("pid_60ms_885",   2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("dt_TOF_60ms_885",2000,tofbd,tofbd+20000,tof, 6000,0,60,decaytime);  
              if(i2>10000 && i2<16000){
                tof2 = tof - partof2*i2; 
                FillHistogram("pid_60ms_885_tofCor",   2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
                FillHistogram("dt_TOF_60ms_885_tofCor",2000,tofbd,tofbd+20000,tof2, 6000,0,60,decaytime);  
              }
            }
          }
          if(e>=900 && e<=906){//BGR_885keV
            if(decaytime<60){
              FillHistogram("pid_60ms_BGR885",   2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("dt_TOF_60ms_BGR885",2000,tofbd,tofbd+20000,tof, 6000,0,60,decaytime);  
              if(i2>10000 && i2<16000){
                tof2 = tof - partof2*i2; 
                FillHistogram("pid_60ms_BGR885_tofCor",   2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
                FillHistogram("dt_TOF_60ms_BGR885_tofCor",2000,tofbd,tofbd+20000,tof2, 6000,0,60,decaytime);  
              }
            }
          }
        }
      }
    }
    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,entries);
      fflush(stdout);
    }


  }
  printf("   on entry %lu / %lu   \n",x,entries);
  SaveHistograms(Form("beta_prompt_op%s.root",runnum.c_str()));

  return;
}



void Histogram::BetaPID(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation1bestT/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  long x=0;
  long entries = beta->GetEntries();
  double dE, dnum, dt, sumE;

  //========== PID cuts for Ne Chain ===========//
  TList *gList = new TList();
  TCutG *cutg;

  //==== PID cuts for Na chain ====//
  double lna = 700;
  double wna = 800;
  double horna = 230;
  double verna = 40;
  for(int i=0;i<10;i++){
    double start[2] = {15700,6260};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("recna%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]+wna);
    cutg->SetPoint(3,start[0],start[1]+wna);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);
  }
  for(int i=0;i<16;i++){
    double start[2] = {13500,5460};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("rec%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]+wna);
    cutg->SetPoint(3,start[0],start[1]+wna);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);
  }

  TCutG *NaCut = new TCutG("NaCut",5);
  NaCut->SetPoint(0,15200,6300);
  NaCut->SetPoint(1,15200,7300);
  NaCut->SetPoint(2,19500,6600);
  NaCut->SetPoint(3,19500,5600);
  NaCut->SetPoint(4,15200,6300);

  TCutG *NeCut = new TCutG("NeCut",5);
  NeCut->SetPoint(0,13500,5400);
  NeCut->SetPoint(1,13500,6400);
  NeCut->SetPoint(2,19100,5500);
  NeCut->SetPoint(3,19100,4500);
  NeCut->SetPoint(4,13500,5400);

  double partof1 = -0.272716;  //tan(TOF = tan*i2)
  double partof2  = -0.257303; //tan(TOF = tan*i2)
  double tof1, tof2;
  int tofbd = 6000;
  for(x=0;x<entries;x++){
    fbeta->Clear();
    beta->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    double tof = fimp.fI2S;
    double i2 = fimp.fI2S_I2N;
    double pin1e = fimp.fPIN1E;

    FillHistogram("pid",2000,tofbd,tofbd+20000,tof, 300,4500,7500,pin1e);
    FillHistogram("i2_tof", 2000,tofbd,tofbd+20000,tof, 4000,0,40000,i2);
    FillHistogram("tof_i2",4000,0,40000,i2,   2000,tofbd,tofbd+20000,tof);
    FillHistogram("pin1_i2", 4000,0,40000,i2,300,4500,7500,pin1e);
    if(i2>10000 && i2<16000){
      FillHistogram("pid_before",2000,tofbd,tofbd+20000,tof, 300,4500,7500,pin1e);
      FillHistogram("i2_before", 2000,tofbd,tofbd+20000,tof, 2000,5000,25000,i2);
      FillHistogram("tof_before",4000,0,40000,i2,   2000,tofbd,tofbd+20000,tof);
      tof2 = tof - partof2*i2; 
      FillHistogram("pid_Pafter",2000,tofbd,tofbd+20000,tof2, 300,4500,7500,pin1e);
      FillHistogram("i2_Pafter", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,i2);
      FillHistogram("tof_Pafter",3000,0,30000,i2,             2000,tofbd,tofbd+20000,tof2);
    }
    for(int d=0;d<fbeta->DecaySize();d++){
      double decaytime = fbeta->fDecay[d].fDecayTime;
      if(decaytime>1000) continue;
      for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
        double e = fbeta->fDecay[d].fGe[m].GetEnergy();
        if(e<10 || e>4000) continue;
        FillHistogram("dt_singles",1000,0,1000,decaytime,4000,0,4000,e);
        if(fimp.fI2S_I2N==32768){
          FillHistogram("dt_singles_i2",1000,0,1000,decaytime,4000,0,4000,e);
        }
        if(fimp.fI2S_I2N==0){
          FillHistogram("dt_singles_i20",1000,0,1000,decaytime,4000,0,4000,e);
        }
        if(e>=882 && e<=888){//885keV
          if(decaytime<100){
            FillHistogram("pid_100ms_885W",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
            FillHistogram("I2_TOF_100ms_885W", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
          }
        }
        if(e>=900 && e<=906){//BGR_885keV
          if(decaytime<100){
            FillHistogram("pid_100ms_BGR885W",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
            FillHistogram("I2_TOF_100ms_BGR885W", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
          }
        }
        if(e>=146 && e<=152){//150keV
          if(decaytime<50){
            FillHistogram("pid_50ms_150W",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
            FillHistogram("I2_TOF_50ms_150W", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
          }
        }
        if(e>=156 && e<=162){//BGR_150keV
          if(decaytime<50){
            FillHistogram("pid_50ms_BGR150W",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
            FillHistogram("I2_TOF_50ms_BGR150W", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
          }
        }
        if(i2>10000 && i2<16000){
          if(e>=882 && e<=888){//885keV
            if(decaytime<100){
              FillHistogram("pid_100ms_885",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_885", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_100ms_885_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_885_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
          if(e>=900 && e<=906){//BGR_885keV
            if(decaytime<100){
              FillHistogram("pid_100ms_BGR885",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_BGR885", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_100ms_BGR885_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_BGR885_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
          if(e>=864 && e<=870){//BGL_885keV
            if(decaytime<100){
              FillHistogram("pid_100ms_BGL885",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_BGL885", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_100ms_BGL885_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_100ms_BGL885_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
          if(e>=146 && e<=152){//150keV
            if(decaytime<50){
              FillHistogram("pid_50ms_150",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_150", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_50ms_150_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_150_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
          if(e>=156 && e<=162){//BGR_150keV
            if(decaytime<50){
              FillHistogram("pid_50ms_BGR150",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_BGR150", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_50ms_BGR150_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_BGR150_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
          if(e>=136 && e<=142){//BGL_150keV
            if(decaytime<50){
              FillHistogram("pid_50ms_BGL150",    2000,tofbd,tofbd+20000,tof, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_BGL150", 2000,tofbd,tofbd+20000,tof, 3000,0,30000,fimp.fI2S_I2N);
              FillHistogram("pid_50ms_BGL150_cor",    2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
              FillHistogram("I2_TOF_50ms_BGL150_cor", 2000,tofbd,tofbd+20000,tof2, 3000,0,30000,fimp.fI2S_I2N);
            }
          }
        }
      }
    }


    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,entries);
      fflush(stdout);
    }


  }
  printf("   on entry %lu / %lu   \n",x,entries);
  SaveHistograms(Form("beta_prompt_op%s.root",runnum.c_str()));


}

void Histogram::BetaSort(){


  //=========== Read cross talk pars ============//


  std::map<int,double[4][4]> xtalmat = ReadMat();    
  printf("Read and write ct paras into xtalmat with %i\n",xtalmat.size());
  //for(int a=0;a<xtalmat.size();a++){
  //  printf("%i \t %i \t ",a,a%4);
  //  for(int b=0;b<4;b++){
  //    for(int c=0;c<4;c++){
  //      printf("%f",xtalmat[a][b][c]);
  //    }
  //    printf("\n");
  //  }
  //}
  //============================================//


  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation1bestT/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)cutf->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");
  TCutG *Na = (TCutG *)cutf1->Get("Na");

  std::map<pixel,Implant> fImpMap;
  std::map<pixel, int> fDecaySize;
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fDecaySize[pix] = -1;
    }
  }

  long x=0;
  long entries = beta->GetEntries();
  double dE, dnum, dt, sumE;


  TCutG *cutg;
  //==== PID cuts====//
  double lna = 300;
  double wna = 1000;
  double horna = 100;
  double verna = 20;
  for(int i=0;i<18;i++){
    double start[2] = {12500,6200};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("recna%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]+wna);
    cutg->SetPoint(3,start[0],start[1]+wna);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);
  }
  for(int i=0;i<22;i++){
    double start[2] = {10770,6140};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("rec%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]-wna);
    cutg->SetPoint(3,start[0],start[1]-wna);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);
  }

  TCutG *NeCut = new TCutG("NeCut1",5);
  NeCut->SetPoint(0,10970,5100);
  NeCut->SetPoint(1,11970,5000);
  NeCut->SetPoint(2,11970,5000+wna);
  NeCut->SetPoint(3,10970,5100+wna);
  NeCut->SetPoint(4,10970,5100);

  TCutG *NeCut0 = new TCutG("NeCut0",5);
  NeCut0->SetPoint(0,11970,5000);
  NeCut0->SetPoint(1,12870,4820);
  NeCut0->SetPoint(2,12870,4820+wna);
  NeCut0->SetPoint(3,11970,5000+wna);
  NeCut0->SetPoint(4,11970,5000);


  double partof1 = -0.272716;  //tan(TOF = tan*i2)
  double partof2  = -0.257303; //tan(TOF = tan*i2)
  double tof1, tof2;
  int tofbd = 6000;
  for(x=0;x<entries;x++){
    fbeta->Clear();
    beta->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    double tof = fimp.fI2S;
    double i2 = fimp.fI2S_I2N;
    double pin1e = fimp.fPIN1E;

    //if(i2>10000 && i2<16000){
    FillHistogram("pid_before",2000,tofbd,tofbd+20000,tof, 300,4500,7500,pin1e);
    FillHistogram("i2_before", 2000,tofbd,tofbd+20000,tof, 2000,5000,25000,i2);
    FillHistogram("tof_before",4000,0,40000,i2,   2000,tofbd,tofbd+20000,tof);
    //tof2 = tof - partof2*i2; 
    tof2 = tof; 
    /*for(int mov=0;mov<18;mov++){
      TCutG *rec = (TCutG *)gList->FindObject(Form("recna%i",mov));
      if(rec->IsInside(tof2, fimp.fPIN1E)){
      FillHistogram(Form("PID_recna%i",mov),2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
      FillHistogram(Form("I2_TOF_recna%i",mov), 2000,tofbd,tofbd+20000,tof2, 2500,5000,30000,fimp.fI2S_I2N);
      }
      }
      for(int d=0;d<fbeta->DecaySize();d++){
      double decaytime = fbeta->fDecay[d].fDecayTime;
      if(decaytime>100) continue;
      std::map<int, Clover> clmap;
      for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
      double gamEm = fbeta->fDecay[d].fGe[m].GetEnergy();
      if(gamEm<10 || gamEm>4000) continue;
      int clnum = (fbeta->fDecay[d].fGe[m].GetNumber()-208)/4;
      clmap[clnum].Add(fbeta->fDecay[d].fGe[m]);
      for(int mov=0;mov<18;mov++){
      TCutG *rec = (TCutG *)gList->FindObject(Form("recna%i",mov));
      if(rec->IsInside(tof2, fimp.fPIN1E)){
      FillHistogram(Form("dt_singles_recna%i",mov), 10000,0,100,decaytime, 1000,0,1000,gamEm);
      }
      }
      }
      }*/
    //}//end of Na implant gate;

    if(Na32->IsInside(tof2, fimp.fPIN1E)){
      FillHistogram(Form("PID_%s",Na32->GetName()),2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
      FillHistogram(Form("I2_TOF_%s",Na32->GetName()), 2000,tofbd,tofbd+20000,tof2, 2500,5000,30000,fimp.fI2S_I2N);
    }
    if(NeCut->IsInside(tof2, fimp.fPIN1E)){
      FillHistogram(Form("PID_%s",NeCut->GetName()),2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
      FillHistogram(Form("I2_TOF_%s",NeCut->GetName()), 2000,tofbd,tofbd+20000,tof2, 2500,5000,30000,fimp.fI2S_I2N);
    }
    if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
      FillHistogram(Form("PID_%s",NeCut0->GetName()),2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
      FillHistogram(Form("I2_TOF_%s",NeCut0->GetName()), 2000,tofbd,tofbd+20000,tof2, 2500,5000,30000,fimp.fI2S_I2N);
    }
    /*for(int mov=0;mov<22;mov++){
      TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
      if(rec->IsInside(tof2, fimp.fPIN1E)){
        FillHistogram(Form("PID_rec%i",mov),2000,tofbd,tofbd+20000,tof2, 300,4500,7500,fimp.fPIN1E);
        FillHistogram(Form("I2_TOF_rec%i",mov), 2000,tofbd,tofbd+20000,tof2, 2500,5000,30000,fimp.fI2S_I2N);
      }
    }*/
    for(int d=0;d<fbeta->DecaySize();d++){
      double decaytime = fbeta->fDecay[d].fDecayTime;
      if(Na32->IsInside(tof2,fimp.fPIN1E)){
        FillHistogram(Form("decaytime_%s",Na32->GetName()),10000,0,1000,decaytime);
      }
      if(decaytime>500) continue;
      std::map<int, Clover> clmap;
      for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
        double gamEm = fbeta->fDecay[d].fGe[m].GetEnergy();
        if(gamEm<10 || gamEm>4000) continue;
        int clnum = (fbeta->fDecay[d].fGe[m].GetNumber()-208)/4;
        clmap[clnum].Add(fbeta->fDecay[d].fGe[m]);
        //for(int mov=0;mov<22;mov++){
        //  TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
        //  if(rec->IsInside(tof2, fimp.fPIN1E)){
        //    FillHistogram(Form("dt_singles_rec%i",mov), 10000,0,100,decaytime, 1000,0,1000,gamEm);
        //  }
        //}
        if(Na32->IsInside(tof2, fimp.fPIN1E)){
          FillHistogram(Form("dt_singles_%s",Na32->GetName()), 500,0,500,decaytime, 4000,0,4000,gamEm);
        }
        for(int n=m+1;n<fbeta->fDecay[d].GeSize();n++){
          double gamEn = fbeta->fDecay[d].fGe[n].GetEnergy();
          if(gamEn<10 || gamEn>4000) continue;
          if(decaytime<50){
            FillHistogram(Form("ggmat50ms_%s",Na32->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
            FillHistogram(Form("ggmat50ms_%s",Na32->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
          }
          if(decaytime<70){
            FillHistogram(Form("ggmat70ms_%s",Na32->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
            FillHistogram(Form("ggmat70ms_%s",Na32->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
          }
        }
        if(NeCut->IsInside(tof2, fimp.fPIN1E)){
          FillHistogram(Form("dt_singles_%s",NeCut->GetName()), 500,0,500,decaytime, 4000,0,4000,gamEm);
        }
        if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
          FillHistogram(Form("dt_singles_%s",NeCut0->GetName()), 500,0,500,decaytime, 4000,0,4000,gamEm);
        }
        for(int n=m+1;n<fbeta->fDecay[d].GeSize();n++){
          double gamEn = fbeta->fDecay[d].fGe[n].GetEnergy();
          if(gamEn<10 || gamEn>4000) continue;
          if(decaytime<=10){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat10ms_%s",NeCut->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat10ms_%s",NeCut->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat10ms_%s",NeCut0->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat10ms_%s",NeCut0->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
          }
          if(decaytime<=30){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat30ms_%s",NeCut->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat30ms_%s",NeCut->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat30ms_%s",NeCut0->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat30ms_%s",NeCut0->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
          }
          if(decaytime<=50){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat50ms_%s",NeCut->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat50ms_%s",NeCut->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat50ms_%s",NeCut0->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat50ms_%s",NeCut0->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
          }
          if(decaytime<=70){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat70ms_%s",NeCut->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat70ms_%s",NeCut->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmat70ms_%s",NeCut0->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmat70ms_%s",NeCut0->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
          }
          if(decaytime>70 && decaytime<100){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmatbg100ms_%s",NeCut->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmatbg100ms_%s",NeCut->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("ggmatbg100ms_%s",NeCut0->GetName()), 4000,0,4000,gamEm, 4000,0,4000,gamEn);
              FillHistogram(Form("ggmatbg100ms_%s",NeCut0->GetName()), 4000,0,4000,gamEn, 4000,0,4000,gamEm);
            }
          }
        }//end ggmat
      }
      // Addback//
      std::map<int, Clover>::iterator it1;
      for(it1=clmap.begin();it1!=clmap.end();it1++){
        it1->second.SetAddE();
        if(Na32->IsInside(tof2, fimp.fPIN1E)){
          FillHistogram(Form("dt_addback_%s",Na32->GetName()), 100,0,100,decaytime, 4000,0,4000,it1->second.AddbackE());
        }
        //if(NeCut->IsInside(tof2, fimp.fPIN1E)){
        //  FillHistogram(Form("dt_addback_%s",NeCut->GetName()), 100,0,100,decaytime, 4000,0,4000,it1->second.AddbackE());
        //}
        //if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
        //  FillHistogram(Form("dt_addback_%s",NeCut0->GetName()), 100,0,100,decaytime, 4000,0,4000,it1->second.AddbackE());
        //}
        std::map<int, Clover>::iterator it2;
        for(it2=next(it1,1);it2!=clmap.end();it2++){
          it2->second.SetAddE();
          if(Na32->IsInside(tof2,fimp.fPIN1E)){
            if(decaytime<=50){
              FillHistogram(Form("adbmat50ms_%s",Na32->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat50ms_%s",Na32->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(decaytime<=70){
              FillHistogram(Form("adbmat70ms_%s",Na32->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat70ms_%s",Na32->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
          if(decaytime<=10){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat10ms_%s",NeCut->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat10ms_%s",NeCut->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat10ms_%s",NeCut0->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat10ms_%s",NeCut0->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
          if(decaytime<=30){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat30ms_%s",NeCut->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat30ms_%s",NeCut->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat30ms_%s",NeCut0->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat30ms_%s",NeCut0->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
          if(decaytime<=50){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat50ms_%s",NeCut->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat50ms_%s",NeCut->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat50ms_%s",NeCut0->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat50ms_%s",NeCut0->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
          if(decaytime<=70){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat70ms_%s",NeCut->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat70ms_%s",NeCut->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmat70ms_%s",NeCut0->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmat70ms_%s",NeCut0->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
          if(decaytime>=90 && decaytime<100){
            if(NeCut->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmatbg10ms_%s",NeCut->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmatbg10ms_%s",NeCut->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
            if(NeCut0->IsInside(tof2, fimp.fPIN1E)){
              FillHistogram(Form("adbmatbg10ms_%s",NeCut0->GetName()), 4000,0,4000,it1->second.AddbackE(), 4000,0,4000,it2->second.AddbackE());
              FillHistogram(Form("adbmatbg10ms_%s",NeCut0->GetName()), 4000,0,4000,it2->second.AddbackE(), 4000,0,4000,it1->second.AddbackE());
            }
          }
        }  
        //FillHistogram("dt_ab", 1000,0,1000,decaytime, 4000,0,4000, it1->second.AddbackE());
        //for(int mov=0;mov<16;mov++){
        //  TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
        //  if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
        //    FillHistogram(Form("dt_ab_rec%i",mov), 1000,0,1000,decaytime, 4000,0,4000, it1->second.AddbackE());
        //  }
        //}
      }///addback end
    }

    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,entries);
      fflush(stdout);
    }


  }
  printf("   on entry %lu / %lu   \n",x,entries);
  SaveHistograms(Form("beta_prompt_op%s.root",runnum.c_str()));


}



//========================== Event File Sort ===========================//
// Get pid and momentum before and after TOF correction//

std::map<int,double[4][4]> Histogram::ReadMat(std::string filename){

  std::map<int, double[4][4]> xtalmat;
  std::ifstream infile;
  std::string line;
  infile.open(filename.c_str());
  int i=0;
  while(getline(infile,line)){
    int det;
    double c0; 
    double c1; 
    double c2; 
    double c3; 

    std::stringstream ss(line);
    ss >> det;
    ss >> c0;
    ss >> c1;
    ss >> c2;
    ss >> c3;

    xtalmat[det][i%4][0] = c0;
    xtalmat[det][i%4][1] = c1;
    xtalmat[det][i%4][2] = c2;
    xtalmat[det][i%4][3] = c3;
    i +=1;
  }

  return xtalmat;
}

void Histogram::EventSort(){
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  //std::map<int,double[4][4]> xtalmat = ReadMat();

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(fevent->Pin1E()>0){
      FillHistogram("pid",2e3,0,2e4,fevent->I2S(),4e3,0,8e3,fevent->Pin1E());
      if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
        FillHistogram("pid_dssdlo",2e3,0,2e4,fevent->I2S(),4e3,0,8e3,fevent->Pin1E());
      }else{ // no dssd lo
        if(fevent->HGFSize()>0 || fevent->HGBSize()>0){
          FillHistogram("pid_dssdhi",2e3,0,2e4,fevent->I2S(),4e3,0,8e3,fevent->Pin1E());
        }else{ // no dssd hi
          FillHistogram("pid_nodssd",2e3,0,2e4,fevent->I2S(),4e3,0,8e3,fevent->Pin1E());
        }
      }
    }

    for(auto &it:fevent->fHits){
      FillHistogram("sum", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
      FillHistogram("sumsize", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
      if(fevent->Pin1E()>0){
        FillHistogram("sum_pin1", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
        FillHistogram("sumsize_pin1", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
        if(fevent->LGFSize()>0 || fevent->LGBSize()>0) {
          FillHistogram("sum_pin1_dssdlo", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
          FillHistogram("sumsize_pin1_dssdlo", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
        }else{ // no dssd lo
          if(fevent->HGFSize()>0 || fevent->HGBSize()>0){
            FillHistogram("sum_pin1_dssdhi", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
            FillHistogram("sumsize_pin1_dssdhi", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
          }else{ // no dssd hi
            FillHistogram("sum_pin1_nodssd", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
            FillHistogram("sumsize_pin1_nodssd", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
          }
        }
      }else{ // no pin1
        if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
          FillHistogram("sum_nopin1_dssdlo", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
          FillHistogram("sumsize_nopin1_nodssdlo", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
        }else{ // no dssd lo
          if(fevent->HGFSize()>0 || fevent->HGBSize()>0){
            FillHistogram("sum_dec", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
            FillHistogram("sumsize_dec", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
          }else{ // no dssd hi
            if(fevent->HPGeSize()>0 || fevent->LaBr().size()>0){
              FillHistogram("sum_gamma", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
              FillHistogram("sumsize_gamma", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
            }else{ // no gamma
              FillHistogram("sum_sssdonly", 20e3,0,40e3,it.GetEnergy(), 300,0,300,it.GetNumber());
              FillHistogram("sumsize_sssdonly", 100,0,100,fevent->Size(), 300,0,300,it.GetNumber());
            }
          }
        }
      }
      for(auto &it1:fevent->fHits){
        FillHistogram("hitpad_pin1",300,0,300,it.GetNumber(), 300,0,300,it1.GetNumber());
        FillHistogram("hitpad_pin1",300,0,300,it1.GetNumber(), 300,0,300,it.GetNumber());
        FillHistogram("hitpad_pin1E",300,0,300,it.GetNumber(), 300,0,300,it1.GetNumber(),4e3,0,16e3,fevent->Pin1E());
        FillHistogram("hitpad_pin1E",300,0,300,it1.GetNumber(), 300,0,300,it.GetNumber(),4e3,0,16e3,fevent->Pin1E());
      }
    }



    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("event_op5us_%s.root",runnum.c_str()));
  return;


}




//=========================== ListSort ===================================//

void Histogram::Process(std::vector<DetHit> vec, std::map<int,double[4][4]> xtalmat){
  int sssd_flag = 0;
  std::vector<DetHit> FL;
  std::vector<DetHit> BL;
  std::vector<DetHit> Ge;
  DetHit pin1_i2s;
  DetHit pin1;
  DetHit i2s_i2n;
  int num = 0;
  int numy = 0;
  int numz = 0;
  for(size_t y=0;y<vec.size();y++){
    FillHistogram("Summary",16e3,0,16e3,vec[y].GetCharge(), 300,0,300,vec[y].GetNumber());
    FillHistogram("Summary_cal",16e3,0,16e3,vec[y].GetEnergy(), 300,0,300,vec[y].GetNumber());
    switch(vec[y].GetNumber()){
      case 40 ... 79:  //FL
        FL.push_back(vec[y]);
        break;

      case 120 ... 159:  //BL
        BL.push_back(vec[y]);
        break;

      case 160 ... 175:  //SSSD
        if(vec[y].GetCharge()>100){
          sssd_flag++;
        }
        break;

      case 177: //pin1_i2s
        pin1_i2s = vec[y];
        break;

      case 180: //i2s_i2n
        i2s_i2n = vec[y];
        break;

      case 181: //pin1
        pin1 = vec[y];
        break;

      case 208 ... 271: //HPGe
        Ge.push_back(vec[y]);
        break;

      default:
        break;
    }
  }
  if(pin1.GetCharge()>100){
    FillHistogram("PID",2e3,0,2e4,pin1_i2s.GetCharge(),10e3,0,20e3,pin1.GetCharge());
    FillHistogram("Pin1E",10e3,0,20e3,pin1.GetCharge());
  }
  std::map<int,Clover> clmap; 
  if(Ge.size()>0){
    for(auto it:Ge){
      num = it.GetNumber()-208;
      if(it.GetEnergy()>10 && it.GetEnergy()<4000){
        int clnum = num/4;
        int xtalnum = num%4;
        if(clnum == 8) {FillHistogram(Form("single_cl8_%i",xtalnum),8e3,0,4e3,it.GetEnergy());}
        if(clnum == 9) {FillHistogram(Form("single_cl9_%i",xtalnum),8e3,0,4e3,it.GetEnergy());} 
        //if(it.GetEnergy()>1450)
        //printf("channum = [%i][%i] \t Charge = %f \t Energy = %f\n\n",clnum,xtalnum,it.GetCharge(),it.GetEnergy());
        clmap[clnum].Add(it);
        FillHistogram("single",8000,0,4000,it.GetEnergy());
        FillHistogram("single_char",8000,0,4000,it.GetCharge());
        FillHistogram(Form("single_%i",clnum),8000,0,4000,it.GetEnergy());
      }  
    }

    std::map<int,Clover>::iterator it;
    for(it=clmap.begin();it!=clmap.end();it++){
      it->second.SetAddE();
      double adde = 0;
      for(size_t m=0;m<it->second.Size();m++){
        adde += it->second.fXtal[m].GetEnergy();
        for(size_t z=0;z<it->second.Size();z++){
          numy = (it->second.fXtal[m].GetNumber()-208) % 4;
          numz = (it->second.fXtal[z].GetNumber()-208) % 4;
          adde += xtalmat[it->first][numy][numz]*it->second.fXtal[z].GetEnergy();
        }
      }
      FillHistogram("Addback",8000,0,4000,it->second.AddbackE());
      FillHistogram(Form("Addback_%i",it->first),8000,0,4000,it->second.AddbackE());
      it->second.SetAddE(adde);
      FillHistogram("Addback_ct",8000,0,4000,it->second.AddbackE());
      FillHistogram(Form("Addback_ct_%i",it->first),8000,0,4000,it->second.AddbackE()); 
    }
  }
}


void Histogram::ListSort(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  int faddress;
  int fnumber;
  double ftimestamp;
  double fcharge;
  gChain->SetBranchAddress("address",   &faddress);
  gChain->SetBranchAddress("number",    &fnumber);
  gChain->SetBranchAddress("timestamp", &ftimestamp);
  gChain->SetBranchAddress("charge",    &fcharge);

  TChannel::ReadDetMapFile();
  std::map<int,double[4][4]> xtalmat = ReadMat();

  double starttime = 0;
  double buildtime = BUILDTIME;
  std::vector<DetHit> vec_hit;

  long n = gChain->GetEntries();
  //n = 1e6;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    DetHit *fhit = new DetHit;
    fhit->SetAddress(faddress);
    fhit->SetNumber(fnumber);
    fhit->SetTimestamp(ftimestamp);
    fhit->SetCharge(fcharge);
    if((fhit->GetTimestamp()-starttime) > buildtime){
      starttime = fhit->GetTimestamp();
      Process(vec_hit, xtalmat);
      vec_hit.clear();
    }
    vec_hit.push_back(*fhit);
    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("list_output%s.root",runnum.c_str()));
  return;

}



//================================ LogX  ================================//
void Histogram::LogX(){
  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  int Nbins = 1000;
  double xlow = 4; // unit: nanosecond
  double xhigh = 1e9; //5e8 ns = 500ms
  double dx = log(xhigh/xlow)/(double)Nbins;
  Double_t edges[Nbins+1];
  edges[0] = 0;
  for(int i=0;i<Nbins+1;i++){
    edges[i+1] = exp(log(xlow) + (double)i*dx);
  }

  long x=0;
  long n = gChain->GetEntries();

  std::vector<double> tofpara;
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;
  if(tofpara.empty()) return;

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newsubblobNa.root"); // subcuts from pid only Na32 after TOF fluctuation.
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  std::vector<TCutG *>cut2;
  TFile *mycut2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/decaytimecuts.root"); // cuts from decaytime and HPGe.Energy; timecuts[0] is "LHE" cut
  TIter keys2 (mycut2->GetListOfKeys());
  while(TKey *key2 = (TKey*)keys2.Next()) {
    cut2.push_back((TCutG*)key2->ReadObj());
  }
  TH1D *hist = new TH1D("Logx_dtime", "LogX_dtime", Nbins,edges);
  TH1D *hist1 = new TH1D("Logx_dtime1", "LogX_dtime1", Nbins,edges);
  TH1D *hist2 = new TH1D("Logx_dtime2", "LogX_dtime2", Nbins,edges);
  TH1D *hist3 = new TH1D("Logx_dtime3", "LogX_dtime3", Nbins,edges);
  TH1D *hist4 = new TH1D("Logx_dtime4", "LogX_dtime4", Nbins,edges);

  TH1D *ht = new TH1D("Logx_decaytime", "LogX_decaytime", Nbins,edges);
  TH1D *ht1 = new TH1D("Logx_decaytime1", "LogX_decaytime1", Nbins,edges);
  TH1D *ht2 = new TH1D("Logx_decaytime2", "LogX_decaytime2", Nbins,edges);
  TH1D *ht3 = new TH1D("Logx_decaytime3", "LogX_decaytime3", Nbins,edges);
  TH1D *ht4 = new TH1D("Logx_decaytime4", "LogX_decaytime4", Nbins,edges);

  TFile *newf = new TFile(Form("logdecay_%s.root",runnum.c_str()),"recreate");

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    double charge = fbeta->fImplant.fI2S;
    charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
    if(fbeta->fImplant.Stopped()){
      int count = 0;
      for(auto &it1: fbeta->fDecay){
        count++;
        for(auto &it2: it1.fGe){
          double dt = it1.GetTimestamp() - it2.GetTimestamp();
          for(size_t m=0;m<cut1.size();m++){
            if(m==1) {
              if(cut1[m]->IsInside(charge ,fbeta->fImplant.fPIN1E)){
                double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
                dtime = fabs(dtime); //unit: ns;
                if(count<2){
                  hist1->Fill(dtime);
                }
                if(count<3){
                  hist2->Fill(dtime);
                }
                if(count<4){
                  hist3->Fill(dtime);
                }
                if(count<5){
                  hist4->Fill(dtime);
                }
                hist->Fill(dtime);

                if(cut2[0]->IsInside(dt, it2.GetEnergy())){
                  if(count<2){
                    ht1->Fill(dtime);
                  }
                  if(count<3){
                    ht2->Fill(dtime);
                  }
                  if(count<4){
                    ht3->Fill(dtime);
                  }
                  if(count<5){
                    ht4->Fill(dtime);
                  }
                  ht->Fill(dtime);
                }


              }
            }          
          }

        }
      }
    }
    if((x%20000)==0){
      printf("on entry %lu / %lu  \r",x,n);
      fflush(stdout);
    }
  }
  printf("      on entry %lu   /   %lu    \n",x,n);
  hist1->Write();
  hist2->Write();
  hist3->Write();
  hist4->Write();
  hist->Write();
  ht->Write();
  ht1->Write();
  ht2->Write();
  ht3->Write();
  ht4->Write();
  newf->Close();
  return;
}


//========================================== ================================//






