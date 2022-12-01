//g++ junk.cxx `root-config --cflags --libs` -ltbb -L../lib  -lHIT -lBCSUtil  -I../include


#include<cstdio>
#include<string>
#include<iostream>

#include <DDASEvent.h>
#include<DetHit.h>
#include<TChannel.h>
#include<util.h>

#include<TChain.h>
#include<TFile.h>

#define BUILDTIME 5000

TFile  *outfile = 0;
TH2D *pid       = 0;

TCutG *PIDGATE = 0;
TCutG *banana1 = 0;
TCutG *banana2 = 0;


void handleEvent(std::vector<DetHit> &vecHit) {
  double pin1E = -1;
  double tof   = -1;
  bool mpin1E = false;
  bool mtof   = false;


  for(size_t i=0;i<vecHit.size();i++) {
    if(vecHit.at(i).GetNumber()==181) {
      if(pin1E<0) 
        pin1E = vecHit.at(i).GetCharge();
      else
       mpin1E = true;
       //std::cout << "FOUND MULTIPLE PIN1s!!!" << std::endl;
    }
    
    if(vecHit.at(i).GetNumber()==177) {
      if(tof<0) 
        tof = vecHit.at(i).GetCharge();
      else
       mtof = true;
       //std::cout << "FOUND MULTIPLE TOFs!!!" << std::endl;
    }
  }
  
  if(pin1E>0 && tof>0) {
    if(PIDGATE->IsInside(tof,pin1E)) {
      if(!pid) {
        outfile->cd();
        pid = new TH2D("pid","pid",8000,0,0,8000,0,0);
      }
      //std::cout << pin1E << std::endl;
      //std::cout << tof << std::endl;
      if(mpin1E) 
        std::cout << "FOUND MULTIPLE PIN1s!!!" << std::endl;
  
      if(mtof) 
        std::cout << "FOUND MULTIPLE TOFs!!!" << std::endl;

      pid->Fill(tof,pin1E);
    }
  }

  // check pid gate
  // fill a tree

  return;
}


int main(int argc, char **argv) {

  TChain *chain = new TChain("dchan");
  if(argc==1) {
    chain->Add("data/run1040-01.root");
  } else {
    for(int i=1;i<argc;i++) 
      chain->Add(argv[i]);
  }  

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/overview_pidcut_tofexd.root");
  PIDGATE = (TCutG *)cutf1->Get("pid");
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/vetogate.root");
  banana1 = (TCutG *)cutf2->Get("banana1");
  banana2 = (TCutG *)cutf2->Get("banana2");





  TChannel::ReadDetMapFile();

  DDASEvent *event  = new DDASEvent;
  ddaschannel *chan = new ddaschannel;

  chain->SetBranchAddress("ddasevent",&event);

  long nentries = chain->GetEntries();
  long x;
  std::string runnum = GetRunNumber(chain->GetCurrentFile()->GetName());
  outfile = new TFile(Form("outfile%s.root",runnum.c_str()),"RECREATE");

  double starttime = 0;

  std::vector<DetHit> vecHit;

  while(x<nentries) {
    chain->GetEntry(x++);

    for(size_t y=0;y<event->GetNEvents();y++) {
      
      chan = event->GetData()[y];
      DetHit hit(chan);


      if((hit.GetTimestamp()-starttime)>BUILDTIME){
        handleEvent(vecHit);
        starttime = hit.GetTimestamp();
        vecHit.clear();
      }
      vecHit.push_back(hit);
    }

    if((x%5000)==0) {
      printf("   on entry %lu / %lu             \r",x,nentries);
      fflush(stdout);
    }
  }
  printf("   on entry %lu / %lu             \n",x,nentries);
  
  outfile->Write();


  return 0;
}
