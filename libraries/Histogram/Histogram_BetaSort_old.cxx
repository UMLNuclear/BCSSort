
#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<map>
#include<string>

#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
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
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_tofcor/correlation1bestT/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)cutf->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");

  std::map<pixel,Implant> fImpMap;
  std::map<pixel, int> fDecaySize;
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fDecaySize[pix] = -1;
    }
  }

  long x=0;
  long n = beta->GetEntries();
  double dE, dnum, dt, sumE;

  //========== 2D Moving Rectangle1 ===========//
  TList *gList = new TList();
  TCutG *cutg;
  double l = 600;
  double w = 700;
  double hor = 85;
  double ver= 90;
  double hor2 = 60;
  double ver2= 240;


  for(int i=0;i<16;i++){
    double start[2] = {11100,11100};
    start[0] += i*hor;
    start[1] += i*ver;
    cutg = new TCutG(Form("rec%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+l,start[1]-300);
    cutg->SetPoint(2,start[0]+l,start[1]+w-300);
    cutg->SetPoint(3,start[0],start[1]+w);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);  
  }
  for(int i=0;i<10;i++){
    double start[2] = {12400,11100};
    start[0] -= i*hor2;
    start[1] += i*ver2;
    cutg = new TCutG(Form("recr%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+l,start[1]-300);
    cutg->SetPoint(2,start[0]+l,start[1]+w-300);
    cutg->SetPoint(3,start[0],start[1]+w);
    cutg->SetPoint(4,start[0],start[1]);
    gList->Add(cutg);
  }
  //cutg = new TCutG("candy1",7);
  //cutg->SetPoint(0,11100,11100);
  //cutg->SetPoint(1,11100,11800);
  //cutg->SetPoint(2,11640,12310);
  //cutg->SetPoint(3,12240,12010);
  //cutg->SetPoint(4,12240,11310);
  //cutg->SetPoint(5,11700,10800);
  //cutg->SetPoint(6,11100,11100);
  //gList->Add(cutg);

  //cutg = new TCutG("candy2",7);
  //cutg->SetPoint(0,11730,11695);
  //cutg->SetPoint(1,11730,12395);
  //cutg->SetPoint(2,12450,13075);
  //cutg->SetPoint(3,13050,12975);
  //cutg->SetPoint(4,13050,12075);
  //cutg->SetPoint(5,12330,11395);
  //cutg->SetPoint(6,11730,11695);
  //gList->Add(cutg);


  for(x=0;x<n;x++){
    beta->Clear();
    beta->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    bool flagimp = true;
    FillHistogram("pid",4000,0,20000,fimp.fI2S, 4000,0,8000,fimp.fPIN1E);
    //FillHistogram("decaysize",100,0,100,fbeta->DecaySize());
    //FillHistogram("decaysize_1000",100,0,100,fbeta->DecaySize(1000));
    if(Na32->IsInside(fimp.fI2S, fimp.fPIN1E)){
      FillHistogram("PID_32Na", 4000,0,20000,fimp.fI2S, 4000,0,8000,fimp.fPIN1E);
      for(auto &it:fbeta->fDecay){
        FillHistogram("decaytime_32Na",100000,0,10000,it.fDecayTime);
        if(it.fDecayTime>100) continue;
        for(auto &it2:it.fGe){
          if(it2.GetEnergy()<10 || it2.GetEnergy()>4000) continue;
          FillHistogram("single_dt_32Na",1000,0,100,it.fDecayTime, 4000,0,4000,it2.GetEnergy());
          for(auto &it3:it.fGe){
            if(it2.GetTimestamp() == it3.GetTimestamp()) continue;
            if(it3.GetEnergy()<10 || it3.GetEnergy()>4000) continue;
            FillHistogram("ggmat_100_32Na",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
            if(it.fDecayTime>70) continue;
            FillHistogram("ggmat_70_32Na",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
          } 
        }
      }
    }

    if(Ne->IsInside(fimp.fI2S, fimp.fPIN1E)){ // Gate: Ne chain in PID
      FillHistogram("I2_TOF_Ne",2000,0,20000,fimp.fI2S, 3000,0,30000,fimp.fI2S_I2N);
      if((fimp.fI2S>11100 && fimp.fI2S<13000) && (fimp.fI2S_I2N>8000 && fimp.fI2S_I2N<16000)){ // Gate on I2position and TOF
        FillHistogram("decaysize_31Ne",100,0,100,fbeta->DecaySize());
        FillHistogram("decaysize_31Ne_1000",100,0,100,fbeta->DecaySize(1000));
        FillHistogram("I2_TOF_Ne31", 190,11100,13000,fimp.fI2S, 300,10800,13800,fimp.fI2S_I2N);
        FillHistogram("PID_Ne31", 250,11000,13500,fimp.fI2S, 200,4700,6700,fimp.fPIN1E);
        FillHistogram("decaysize_31Ne_50", 100,0,100,fbeta->DecaySize(50));
        for(int mov=1;mov<16;mov++){
          TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
          if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
            FillHistogram(Form("I2_TOF_rec%i",mov), 190,11100,13000,fimp.fI2S, 300,10800,13800,fimp.fI2S_I2N);
            FillHistogram(Form("PID_rec%i",mov), 250,11000,13500,fimp.fI2S, 200,4700,6700,fimp.fPIN1E);
          }
        }
        //for(int mov=0;mov<10;mov++){
        //  TCutG *rec = (TCutG *)gList->FindObject(Form("recr%i",mov));
        //  if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
        //    FillHistogram(Form("I2_TOF_recr%i",mov), 190,11100,13000,fimp.fI2S, 300,10800,13800,fimp.fI2S_I2N);
        //    FillHistogram(Form("PID_recr%i",mov), 250,11000,13500,fimp.fI2S, 200,4700,6700,fimp.fPIN1E);
        //  }
        //}
        for(auto &it:fbeta->fDecay){
          if(it.fDecayTime>1000) continue;
          //FillHistogram("GeSize",64,0,64,it.GeSize());
          std::map<int, Clover> clmap;
          for(auto &it2:it.fGe){
            if(it2.GetEnergy()<10 || it2.GetEnergy()>4000) continue;
            int clnum = (it2.GetEnergy()-208)/4.;
            clmap[clnum].Add(it2);
            FillHistogram("dt_singles", 1000,0,1000,it.fDecayTime, 4000,0,4000, it2.GetEnergy());
            for(int mov=1;mov<16;mov++){
              TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
              if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
                FillHistogram(Form("dt_singles_rec%i",mov), 1000,0,1000,it.fDecayTime, 4000,0,4000, it2.GetEnergy());
              }
            }
            //for(int mov=0;mov<10;mov++){
            //  TCutG *rec = (TCutG *)gList->FindObject(Form("recr%i",mov));
            //  if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
            //    FillHistogram(Form("dt_singles_recr%i",mov), 2000,0,1000,it.fDecayTime, 4000,0,4000, it2.GetEnergy());
            //  }
            //}
            for(auto &it3:it.fGe){
              if(it3.GetEnergy()<10 || it3.GetEnergy()>4000) continue;
              if(it3.GetTimestamp() == it2.GetTimestamp()) continue;
              if(it.fDecayTime<300){
                FillHistogram("ggmat_300",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              if(it.fDecayTime<250){
                FillHistogram("ggmat_250",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              if(it.fDecayTime<200){
                FillHistogram("ggmat_200",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              if(it.fDecayTime<150){
                FillHistogram("ggmat_150",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              if(it.fDecayTime<100){
                FillHistogram("ggmat_100",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              if(it.fDecayTime<50){
                FillHistogram("ggmat_50",4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              }
              for(int mov=1;mov<16;mov++){
                TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
                if(rec->IsInside(fimp.fI2S,fimp.fI2S_I2N)){
                  if(it.fDecayTime<300){
                    FillHistogram(Form("ggmat300_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                  if(it.fDecayTime<250){
                    FillHistogram(Form("ggmat250_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                  if(it.fDecayTime<200){
                    FillHistogram(Form("ggmat200_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                  if(it.fDecayTime<150){
                    FillHistogram(Form("ggmat150_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                  if(it.fDecayTime<100){
                    FillHistogram(Form("ggmat100_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                  if(it.fDecayTime<50){
                    FillHistogram(Form("ggmat50_rec%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                  }
                }
              }            
              //for(int mov=0;mov<10;mov++){
              //  TCutG *rec = (TCutG *)gList->FindObject(Form("recr%i",mov));
              //  if(rec->IsInside(fimp.fI2S,fimp.fI2S_I2N)){
              //    if(it.fDecayTime<50){
              //      FillHistogram(Form("ggmat50_recr%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              //    }
              //    if(it.fDecayTime<30){
              //      FillHistogram(Form("ggmat30_recr%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              //    }
              //    if(it.fDecayTime<15){
              //      FillHistogram(Form("ggmat15_recr%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              //    }
              //    if(it.fDecayTime<10){
              //      FillHistogram(Form("ggmat10_recr%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              //    }
              //    if(it.fDecayTime<5){
              //      FillHistogram(Form("ggmat5_recr%i",mov),4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
              //    }
              //  }
              //}            
            }
          }
          //== Addback ==//
          std::map<int, Clover>::iterator it4;
          for(it4=clmap.begin();it4!=clmap.end();it4++){
            it4->second.SetAddE();
            FillHistogram("dt_ab", 1000,0,1000,it.fDecayTime, 4000,0,4000, it4->second.AddbackE());
            for(int mov=1;mov<16;mov++){
              TCutG *rec = (TCutG *)gList->FindObject(Form("rec%i",mov));
              if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
                FillHistogram(Form("dt_ab_rec%i",mov), 1000,0,1000,it.fDecayTime, 4000,0,4000, it4->second.AddbackE());
              }
            }
            //for(int mov=0;mov<10;mov++){
            //  TCutG *rec = (TCutG *)gList->FindObject(Form("recr%i",mov));
            //  if(rec->IsInside(fimp.fI2S, fimp.fI2S_I2N)){
            //    FillHistogram(Form("dt_ab_recr%i",mov), 2000,0,1000,it.fDecayTime, 4000,0,4000, it4->second.AddbackE());
            //  }
            //}
          }
        }
      }
    }


    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }


  }
  printf("   on entry %lu / %lu   \n",x,n);
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






