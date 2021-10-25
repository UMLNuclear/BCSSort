#include<map>


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



void eventsort(){


  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  //std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);



  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    for(auto &it1:fevent->LGF()){
      FillHistogram("Summary_FL_befor",8000,0,8000,it1.GetCharge(), 40,40,80,it1.GetNumber());
      FillHistogram("Summary_FL_after",4000,0,8000,it1.GetEnergy(),  40,40,80,it1.GetNumber());
    }
    for(auto &it2:fevent->LGB()){
      FillHistogram("Summary_BL_befor",8000,0,8000,it2.GetCharge(), 40,120,160,it2.GetNumber());
      FillHistogram("Summary_BL_after",4000,0,8000,it2.GetEnergy(),  40,120,160,it2.GetNumber());
    }
    for(auto &it3:fevent->fHits){
      FillHistogram("Summary_befor",8000,0,34000,it3.GetCharge(), 300,0,300,it3.GetNumber());
      FillHistogram("Summary_after",4000,0,8000,it3.GetEnergy(),  300,0,300,it3.GetNumber());
    }
    FillHistogram("pid",2000,0,20000,fevent->I2S(), 4000,0,33000,fevent->Pin1E());
    FillHistogram("PIN1_2",4000,0,32000,fevent->Pin1E(), 4000,0,33000,fevent->Pin2E());
    if(fevent->SSSDSize()==0 && fevent->LGFSize()>0 && fevent->LGBSize()>0){
      FillHistogram("pid_nosssd",2000,0,20000,fevent->I2S(), 4000,0,34000,fevent->Pin1E());
      FillHistogram("PIN1_2_nosssd",4000,0,34000,fevent->Pin1E(), 4000,0,34000,fevent->Pin2E());
    }
    
    

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("eventoutput%s.root",runnum.c_str()));
  return;
}
