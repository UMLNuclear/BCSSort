
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


std::map<pixel, double> fChargeMap; //pix = <FL#, BL#>, DetHit = pin1
std::map<pixel, double> fTimeMap; //pix = <FL#, BL#>, DetHit = pin1


std::pair<int, int> GetPixel(std::vector<DetHit> LGF, std::vector<DetHit> LGB){
  int indexf = (-LGF[0].GetNumber()+80)-1;
  int indexb = (-LGB[0].GetNumber()+160)-1;
  double tempE = LGF[0].GetEnergy();
  for(int z=1;z<LGF.size();z++){
    if(tempE<LGF[z].GetEnergy()){
      tempE = LGF[z].GetEnergy();
      indexf = (-LGF[z].GetNumber()+80)-1;   // from 0~39: FL1->0; FL40->39
    }
  }
  tempE = LGB[0].GetEnergy();
  for(int z=1;z<LGB.size();z++){
    if(tempE<LGB[z].GetEnergy()){
      tempE = LGB[z].GetEnergy();
      indexb = (-LGB[z].GetNumber()+160)-1;
    }
  }
  return std::make_pair(indexf,indexb);
}




void eventsort(){

  //std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());

  std::vector<TCutG* >pidcuts;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/pidcut_48.root");
  TIter keys1 (mycut1->GetListOfKeys());  
  while(TKey *key1 = (TKey*)keys1.Next()){
    pidcuts.push_back((TCutG*)key1->ReadObj());
  }


  for(int m=0;m<40;m++){
    for(int n=n;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fChargeMap[pix] = -1;
      fTimeMap[pix] = -1;
    }
  }


  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long nentries = gChain->GetEntries();
  nentries = 1e5;
  long x = 0;
  int pixf, pixb;
  double bmaxE, fmaxE;
  size_t m,n,c;

  for(x=0;x<nentries;x++){
    gChain->GetEntry(x);
    

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,nentries);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,nentries);
  //SaveHistograms(Form("eventdo%s.root",runnum.c_str()));
  SaveHistograms("hitpad_eventdo0048.root");
  return;
}
