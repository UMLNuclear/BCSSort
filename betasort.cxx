#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<map>


#include<TFile.h>
#include<TCutG.h>
#include<TH2D.h>
#include<TChain.h>
#include<TKey.h>
#include<TList.h>

#include"Implant.h"
#include"util.h"
#include"TChannel.h"


using namespace std;




int main(int argc, char** argv){

  TChain *chain = new TChain("beta");
  for(int i=1;i<argc;i++){
    chain->Add(argv[i]);
    printf("added %s to chian. \n", argv[i]);
  }

  Beta *b = 0;
  chain->SetBranchAddress("beta",&b);
  TChannel::ReadDetMapFile();

  long x=0;
  long n = chain->GetEntries();

  vector<TCutG*> cuts;
  TFile *mycut = TFile::Open("~zhu/notebooks/Ne31/Time_Correlation/mynewcuts.root");
  TIter keys (mycut->GetListOfKeys());
  while(TKey *key = (TKey*)keys.Next()) {
    cuts.push_back((TCutG*)key->ReadObj());
  }

  for(x;x<n;x++){
    chain->GetEntry(x);

    for(long y=0;y<cuts.size();y++){
      if(Inside(cuts[y])){
        double indexf = -(b->fImplant.GetPixel().first+1)-80;
        double indexb = -(b->fImplant.GetPixel().second+1)-160;
        double FLT = b->fImplant.fDSSDFront[indexf].GetTimestamp;
        double BLT = b->fImplant.fDSSDFront[indexb].GetTimestamp;
        double dt = BLT - FLT        FillHistogram("timedif_FBLow",cuts[y],200,-1000,1000,dt);
        FillHistogram("timedif_FBLow",cuts[y],200,-1000,1000,dt);
        FillHistogram("timedif_FBLow_FLE",cuts[y],200,-1000,1000,dt, 8000,0,4000,b->fImplant.fDSSDFront[indexf].GetEnergy());
      }
    }
    if((x%5000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }

  }

  printf("   on entry %lu / %lu   \n",x,n);
  SaveHistograms("gamma.root");

  return 0;

}













