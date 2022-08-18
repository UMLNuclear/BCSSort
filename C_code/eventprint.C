
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




void eventprint(){

  TChannel::ReadDetMapFile();

  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/new_good_imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); 
  veccut.push_back(pin1E_dt);

  ofstream ofile;
  ofile.open("entries_check.txt");
  ofile << "run" << "\t"
    << "imp_entries" << "\t"
    << "eve+2gates" << "\t"
    << "imp+fbl" << "\t"
    << "imp+fbl+2gates" << "\t"
    << "pin1e" << "\t"
    << "tof" << "\t"
    << "pin1t" << "\t"
    << "beta" <<std::endl;


  for(int runnum=1023;runnum<1140;runnum++){
    TChain *cheve = new TChain("event");
    BCSEvent *fevent = new BCSEvent;
    cheve->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/event/event%i*.root",runnum));
    cheve->SetBranchAddress("BCSEvent", &fevent);

    TChain *chimp = new TChain("implant");
    Implant *fimp =new Implant;
    chimp->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/implant/imp_good_TOFexd/implant%i*.root",runnum));
    chimp->SetBranchAddress("Implant", &fimp);

    TChain *beta = new TChain("beta");
    Beta *fbeta = new Beta;
    beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation0bestT/beta%i*.root",runnum));
    beta->SetBranchAddress("Beta", &fbeta);

    long neve = cheve->GetEntries();
    long nimp = chimp->GetEntries();
    long x = 0;
    if(nimp<1) continue;

    double dt, dtt,pin1e,pin1t,tof;
    int c1=0;
    int c2=0;
    int c3=0;
    int c4=0;
    int c5=0;
    int c6=0;

    for(x=0;x<neve;x++){
      fevent->Clear();
      cheve->GetEntry(x);
      if(fevent->Pin1E()>0){
        if(fevent->DSSDloT()>0){
          dt = fevent->DSSDloT() - fevent->Pin1T();
          dt = dt/1000.;
          if(veccut[0]->IsInside(dt, fevent->I2S()) && veccut[1]->IsInside(dt, fevent->Pin1E())){
            c1++; // should be equal to nimp  
            if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
              int c22=0;
              int c33=0;
              int c44=0;
              pin1e = -1;
              tof = -1;
              for(auto &it:fevent->fHits){
                if(it.GetNumber()==181){
                  pin1e = it.GetEnergy();
                  pin1t = it.GetTimestamp();
                }
                if(it.GetNumber()==177){
                  tof = it.GetEnergy();
                }
              }
              dtt = fevent->DSSDloT() - pin1t;
              dtt = dtt/1000.;
              if(pin1e>0 && (!(veccut[1]->IsInside(dt, pin1e)))){
                c22++;
              }
              if(tof>0 && (!(veccut[0]->IsInside(dt, tof)))){
                c33++;
              }
              if((!(veccut[0]->IsInside(dtt,fevent->I2S()))) || (!(veccut[1]->IsInside(dt, fevent->Pin1E())))){
                c33++;
              }
              if(c22>0){c2++;}
              if(c33>0){c3++;}
              if(c44>0){c6++;}
            }
          }
        }
      }
    }

    for(x=0;x<nimp;x++){
      fimp->Clear();
      chimp->GetEntry(x);
      if(fimp->FrontSize()>0 && fimp->BackSize()>0){
        c4++;
        dt = fimp->DSSDloT() - fimp->fPIN1T;
        dt = dt/1000.;
        if(veccut[0]->IsInside(dt, fimp->fI2S) && veccut[1]->IsInside(dt, fimp->fPIN1E)){
          c5++; 
        } 
      }
    }

    ofile << runnum << "\t"
      << nimp << "\t"
      << c1 << "\t"
      << c4 << "\t"
      << c5 << "\t"
      << c2 << "\t"
      << c3 << "\t"
      << c6 << "\t"
      << beta->GetEntries()<<std::endl;

  }

  ofile.close();
  return;

}
