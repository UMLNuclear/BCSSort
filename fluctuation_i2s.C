#include<cstdio>
#include<string>

#include<TFile.h>
#include<TChain.h>

#include<Implant.h>
#include<DetHit.h>
#include<util.h>
#include<TChannel.h>





void fluctuation_i2s(){
 
  string outfile = "tof_"; 
  std::string num = GetRunNumber(_file0->GetName());
  outfile.append(num);
  outfile.append(".root");

  TFile *f = TFile::Open(Form("root_file/tof_histogram/tof%s.root",num.c_str()));
  TH2D *hist = (TH2D *)f->Get("TOF");
  if(hist==NULL) {
    printf("No histogram\n");
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
  //n = 1e5;
  double first_time = -1;
  double runtime;
  double tof,ctof;

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
      runtime = (fevent->Pin1T() - first_time)/1e9;
      tof = fevent->I2S();
      ctof = tof + (offset-sp3->Eval(runtime));
      FillHistogram("ctof", 3800,0,3800,runtime, 1600,0,32000,ctof);
      FillHistogram("tof", 3800,0,3800,runtime, 1600,0,32000,tof);
    }
    if((x%50000)==0){
      printf("on enry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("on enry %lu / %lu \n",x,n);
  SaveHistograms("ctof.root");





}
