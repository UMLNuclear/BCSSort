#include<TFile.h>
#include<TCutG.h>
#include<TH2D.h>
#include<TChain.h>
#include<TKey.h>
#include<TList.h>
#include<TCutG.h>

#include<Implant.h>
#include<util.h>
#include<TChannel.h>
#include<DetHit.h>
#include<TOFCorrection.h>
#include<globals.h>


void eventdt(){

  int runnum = 1040;
  TChain *cheve = new TChain("event");
  BCSEvent *fevent = new BCSEvent;
  cheve->Add(Form("data/5us_tofcor/sample/event%i*",runnum));
  cheve->SetBranchAddress("BCSEvent", &fevent);

  TChannel::ReadDetMapFile();

  long n = cheve->GetEntries();
  //n = 1e5;
  long x = 0;
  double dt;

  for(x=0;x<n;x++){
    cheve->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->DSSDloT()>0){
        dt = fevent->DSSDloT() - fevent->Pin1T();
        dt = dt/1000.;
        FillHistogram("Pin1E_dt", 1000,-5,5,dt,10e2,0,10e3,fevent->Pin1E());
        FillHistogram("TOF_dt", 1000,-5,5,dt,20e2,0,20e3,fevent->I2S());
      }
    }else{
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        dt = fevent->HGFMax().GetTimestamp() - fevent->HGBMax().GetTimestamp();
        dt = dt/1000.;
        FillHistogram("dec_dt_hiFandB", 1000,-5,5,dt);
      }
    }

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("event_op5us_dt_%i_tofcor.root",runnum));
  return;

}


