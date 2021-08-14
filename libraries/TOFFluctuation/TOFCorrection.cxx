#include <cstdio>
#include <string>

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

  TFile *corf = TFile::Open(Form("tof%s.root",num.c_str()));
  TH2D *hist = (TH2D *)corf->Get(Form("tof%s",num.c_str()));
  if(hist==NULL){ //"NULL" = exist;
    printf("histogram does not exist\n");
    return;
  }
  if((hist->Integral())<10){
    printf("histogram is empty\n");
    return;
  }
  TProfile *pfx = hist->ProfileX();
  TSpline3 *sp3 = new TSpline3(pfx);
  TF1 px("px","pol1");
  pfx->Fit(&px);
  double offset = px.Eval(0);

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long x = 0;
  long n = gChain->GetEntries();
  double first_time = -1;
  double runtime;
  double ctof,tof;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
        runtime = (fevent->Pin1T()-first_time)/1e9;
        tof = fevent->I2S();
        ctof = tof + (offset-sp3->Eval(runtime));
        FillHistogram(Form("ctof%s", num.c_str()), 3600,0,3600,runtime, 3000,0,30000,ctof);
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("  on entry %lu / %lu  \n",x,n);
  SaveHistograms(Form("ctof%s.root",num.c_str()));
  
}









