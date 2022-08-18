
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

#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<map>




void eventhitpad(){

  TChain *cheve = new TChain("event");
  BCSEvent *fevent = new BCSEvent;
  cheve->Add("data/5us_tofcor/sample/event104*");
  cheve->SetBranchAddress("BCSEvent", &fevent);

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pidcut104.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }
  
  TChannel::ReadDetMapFile();

  long n = cheve->GetEntries();
  //n = 1e5;
  long x = 0;

  for(x=0;x<n;x++){
    cheve->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
        int lgb = fevent->LGBMax().GetNumber();
        lgb = -lgb + 159;
        int lgf = fevent->LGFMax().GetNumber();
        lgf = -lgf + 79;
        FillHistogram("hitpad_pin1", 40,0,40,lgb, 40,0,40,lgf);
        if(cut1[0]->IsInside(fevent->I2S(), fevent->Pin1E())){
          FillHistogram("hitpad_Na", 40,0,40,lgb, 40,0,40,lgf);
        }
        if(cut1[1]->IsInside(fevent->I2S(), fevent->Pin1E())){
          FillHistogram("hitpad_Ne", 40,0,40,lgb, 40,0,40,lgf);
        }
      }
    }

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms("event5us_hitpad_104_tofcor.root");
  return;

}
