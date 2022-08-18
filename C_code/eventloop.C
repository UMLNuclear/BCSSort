
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
  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/new_good_imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); // TOF vs timedif for 44S;
  veccut.push_back(pin1E_dt); // Pin1E vs timedif for S;

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne"); 

  long n = cheve->GetEntries();
  long x = 0;
  double dt;

  long c1 = 0;


  for(x=0;x<n;x++){
    cheve->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->DSSDloT()>0){
        dt = fevent->DSSDloT() - fevent->Pin1T();
        dt = dt/1000.;
        if(veccut[0]->IsInside(dt, fevent->I2S()) && veccut[1]->IsInside(dt, fevent->Pin1E())){
          c1++;
        }
      }
    }  

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  std::cout<< "good imp = "<<c1<<std::endl;

  return;

}
