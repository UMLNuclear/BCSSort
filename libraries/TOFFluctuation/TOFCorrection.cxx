#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSpline.h>
#include <TChain.h>

#include <util.h>
#include <TChannel.h>
#include <TOFCorrection.h>
#include <DetHit.h>
#include <BCSint.h>



//TChain *gChain = new TChain("event");



TOFCorrection *TOFCorrection::fTOFCorrection = 0;

TOFCorrection *TOFCorrection::Get(){
  if(fTOFCorrection){
    fTOFCorrection = new TOFCorrection;
  }
  return fTOFCorrection;
}

TOFCorrection::TOFCorrection(){}

TOFCorrection::~TOFCorrection(){}

std::vector<double> TOFCorrection::ReadFile(int num, std::string filename){

  std::vector<double> tofpara; 
  std::ifstream infile;
  std::string line;
  infile.open(filename.c_str());
  
  while(getline(infile,line)){
    if(line[0]=='#') continue;
    int runnum;
    double p0;
    double p1;
    double p2;

    std::stringstream ss(line);
    ss >> runnum;
    ss >> p0;
    ss >> p1;
    ss >> p2;

    if(runnum == num){
      tofpara.push_back(p0);
      tofpara.push_back(p1);
      tofpara.push_back(p2);
      return tofpara;
    }
  }
  return tofpara;
}


void TOFCorrection::Fluctuation(){
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  std::string num = runnum.substr(0,4);

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long x = 0;
  long n = gChain->GetEntries();
  double first_time = -1;
  double runtime;
  double tof;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
      runtime = (fevent->Pin1T()-first_time)/1e9;
      tof = fevent->I2S();
      FillHistogram(Form("tof%s",num.c_str()), 3800,0,3800,runtime, 1600,0,32000,tof);
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("  on entry %lu / %lu  \n",x,n);
  SaveHistograms(Form("tof%s.root", num.c_str()));
}

void TOFCorrection::Correct(){

  
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  std::string num = runnum.substr(0,4);


  TFile *corf = TFile::Open(Form("/home/zhu/packages/BCSSort/root_file/tof_correction/tof/tof%s.root",num.c_str()));
  TH2D *hist = (TH2D *)corf->Get(Form("tof%s",num.c_str()));
  if(hist==NULL){ //"NULL" = exist;
    printf("histogram does not exist\n");
    return;
  }
  if((hist->Integral())<10){
    printf("histogram is empty\n");
    return;
  }
  std::vector<double> tofpara = ReadFile(std::stoi(num));
  TProfile *pfx = hist->ProfileX();
  TSpline3 *sp3 = new TSpline3(pfx);
  TF1 px("px","pol1");
  pfx->Fit(&px);
  double offset = px.Eval(0);

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  TH2D *pid = new TH2D("pid_ctof","pid_ctof",2e3,0,2e4, 8e3,0,8e3);
  TH2D *pid2 = new TH2D("pid_newctof","pid_newctof",2e3,0,2e4, 8e3,0,8e3);
  long x = 0;
  long n = gChain->GetEntries();
  double first_time = -1;
  double runtime;
  double ctof2,ctof,tof;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
        runtime = (fevent->Pin1T()-first_time)/1e9;
        tof = fevent->I2S();
        ctof = tof + (offset-sp3->Eval(runtime));
        ctof2 = tofpara[0]+tofpara[1]*ctof+tofpara[2]*ctof*ctof;
        FillHistogram(Form("ctof%s", num.c_str()), 3600,0,3600,runtime, 3000,0,30000,ctof);
        FillHistogram(Form("newctof%s",num.c_str()), 3600,0,3600,runtime, 3000,0,30000,ctof2);
        if(fevent->SSSDSize()==0){
          pid->Fill(ctof, fevent->Pin1E());      
          pid2->Fill(ctof2, fevent->Pin1E());
        }      
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("  on entry %lu / %lu  \n",x,n);
  SaveHistograms(Form("newctof%s.root",num.c_str()));
  TFile *newf = new TFile(Form("PID_TOF%s.root", num.c_str()),"recreate");
  pid->Write(); 
  pid2->Write(); 
  newf->Close();
}









