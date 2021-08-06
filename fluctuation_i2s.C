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
  outfile.append(GetRunNumber(_file0->GetName()));
  outfile.append(".root");

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long x = 0;
  long n = gChain->GetEntries();

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(fevent->I2S()>0){
      for(auto &it1 : fevent->fHits){
        //printf("I2S = %f\n", fevent->I2S());
        FillHistogram("TOF", 3600,0,3600,it1.GetTimestamp()/1e9, 3000,0,30000,fevent->I2S());
      }
    }
    if((x%50000)==0){
      printf("on enry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("on enry %lu / %lu \n",x,n);
  SaveHistograms(outfile);





}
