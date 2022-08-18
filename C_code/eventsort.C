
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






void eventsort(){

  int runnum = 1042;
  TChain *cheve = new TChain("event");
  BCSEvent *fevent = new BCSEvent;
  cheve->Add(Form("data/5us_tofcor/event/event%i*",runnum));
  cheve->SetBranchAddress("BCSEvent", &fevent);

  TChannel::ReadDetMapFile();
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)mycut1->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");


  long n = cheve->GetEntries();
  long x = 0;
  double dt;

  for(x=0;x<n;x++){
    cheve->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->LGFSize()>1){
        DetHit temp = fevent->LGFMax();
        for(int i=0;i<fevent->LGFSize();i++){
          if(temp.GetTimestamp()==fevent->LGF().at(i).GetTimestamp()) continue;
          if(fabs(temp.GetEnergy()-fevent->LGF().at(i).GetEnergy())<2){
            dt = fabs(temp.GetTimestamp()-fevent->LGF().at(i).GetTimestamp());
            dt=dt/1000.;//unit: us
            FillHistogram("LGF_dt_maxE",1000,0,10,dt);
          }
        }
      }
      if(fevent->LGBSize()>1){
        DetHit temp = fevent->LGBMax();
        for(int i=0;i<fevent->LGBSize();i++){
          if(temp.GetTimestamp()==fevent->LGB().at(i).GetTimestamp()) continue;
          if(fabs(temp.GetEnergy()-fevent->LGB().at(i).GetEnergy())<2){
            dt = fabs(temp.GetTimestamp()-fevent->LGB().at(i).GetTimestamp());
            dt=dt/1000.;//unit: us
            FillHistogram("LGB_dt_maxE",1000,0,10,dt);
          }
        }
      }
    }


    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("eventop_%i.root",runnum));
  return;

}
