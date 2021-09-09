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

#include<Implant.h>
#include<util.h>
#include<TChannel.h>


using namespace std;




void betasort (){

  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  long x=0;
  long n = gChain->GetEntries();
  cout<<n<<endl;

  vector<TCutG*> pidcuts;
  vector<TCutG*> timecuts;
  vector<TCutG*> imptimecuts;
  vector<TCutG*> pidnacuts;
  vector<TCutG*> pidnasubcuts;
  vector<TCutG*> pidnecuts;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/blobcuts.root"); //cuts from pid
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    pidcuts.push_back((TCutG*)key1->ReadObj());
  }
  TFile *mycut2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/decaytimecuts.root"); // cuts from decaytime and HPGe.Energy
  TIter keys2 (mycut2->GetListOfKeys());
  while(TKey *key2 = (TKey*)keys2.Next()) {
    timecuts.push_back((TCutG*)key2->ReadObj());
  }
  TFile *mycut3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/implanttimecuts.root"); // cuts from implanttime and HPGe.Energy
  TIter keys3 (mycut3->GetListOfKeys());
  while(TKey *key3 = (TKey*)keys3.Next()) {
    imptimecuts.push_back((TCutG*)key3->ReadObj());
  }
  TFile *mycut4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/PID_Na32_cuts.root"); // cuts from pid only Mg32 
  TIter keys4 (mycut4->GetListOfKeys());
  while(TKey *key4 = (TKey*)keys4.Next()) {
    pidnacuts.push_back((TCutG*)key4->ReadObj());
  }
  TFile *mycut5 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/PID_Na32_subcuts.root"); // cuts from pid only Mg32 
  TIter keys5 (mycut5->GetListOfKeys());
  while(TKey *key5 = (TKey*)keys5.Next()) {
    pidnasubcuts.push_back((TCutG*)key5->ReadObj());
  }
  TFile *mycut6 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/PID_Ne30_cuts.root"); // cuts from pid only Mg32 
  TIter keys6 (mycut6->GetListOfKeys());
  while(TKey *key6 = (TKey*)keys6.Next()) {
    pidnecuts.push_back((TCutG*)key6->ReadObj());
  }

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(fbeta->fImplant.Stopped()){
      int pidflag = 0;
      int pidflag2 = 0;
      for(auto &it1 :fbeta->fDecay){
        if(it1.fImplantTime<50 && !pidflag){
          FillHistogram("pid_50",2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
          pidflag++;
        }
        if(pidflag2>0) continue;
        FillHistogram("pid",2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
        pidflag2++;
      }
    }
    //========= blobcuts from pid============//
    for(long y=0;y<pidcuts.size();y++){ // PID blob
      if(fbeta->fImplant.Stopped()){
        int pidflag = 0;
        if(fbeta->fImplant.Inside(pidcuts[y])){
          double indexf = -(fbeta->fImplant.GetPixel().first+1)-80;
          double indexb = -(fbeta->fImplant.GetPixel().second+1)-160;
          double FLT = fbeta->fImplant.fDSSDFront[indexf].GetTimestamp();
          double BLT = fbeta->fImplant.fDSSDFront[indexb].GetTimestamp();
          double dt = BLT - FLT;
          //FillHistogram("timedif_FBLow",pidcuts[y],200,-1000,1000,dt);
          //FillHistogram("timedif_FBLow_FLE",pidcuts[y],200,-1000,1000,dt, 8000,0,4000,fbeta->fImplant.fDSSDFront[indexf].GetEnergy());
          for(auto &it1: fbeta->fDecay){
            std::map<int, Clover> clmap;
            for(auto &it2 : it1.fGe){
              int clnum = (it2.GetNumber()-208)/4;
              clmap[clnum].Add(it2);
              //double dt_fhge = it1.GetTimestamp() - it2.GetTimestamp(); 
            }
            //========== decaytime cuts =========//
            for(long z=0;z<timecuts.size();z++){ //LHE => timecuts[0]
              std::map<int,Clover>::iterator it;
              for(it=clmap.begin();it!=clmap.end();it++){
                it->second.SetAddE();
                int cflag = 0; // fill addback one time
                for(auto &it3 : it->second.fXtal){
                  double dt_fhge = it1.GetTimestamp() - it3.GetTimestamp(); // FH.t - fGe.t
                  if(timecuts[z]->IsInside(dt_fhge,it3.GetEnergy())){  // LHE => timecuts[0]
                    FillHistogram("summary_single", pidcuts[y], timecuts[z],8000,0,4000,it3.GetEnergy(), 300,0,300,it3.GetNumber());
                    FillHistogram("single", pidcuts[y], timecuts[z],8000,0,4000,it3.GetEnergy());
                    FillHistogram("single_implanttime", pidcuts[y], timecuts[z], 1000,0,500,it1.fImplantTime,8000,0,4000,it3.GetEnergy());
                    for(int m=0;m<imptimecuts.size();m++){ // imptimecuts[0]=>Mg32 peak; imptimecuts[1]=>Na30 peak;
                      if(!pidflag && imptimecuts[m]->IsInside(it1.fImplantTime, it3.GetEnergy())){
                        FillHistogram("pid_nosssd",imptimecuts[m],2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
                        pidflag++;
                      }
                      if(m==0 &&  imptimecuts[m]->IsInside(it1.fImplantTime, it3.GetEnergy())){ // m==0 ==> Na32
                        for(int n=0;n<pidnacuts.size();n++){
                          if(fbeta->fImplant.Inside(pidnacuts[n])){ //blob1 + blob2 from PID_only_Na32 => blob1 is real Mg32; blob2 is compton scattering
                            FillHistogram("gamma_Na32",pidnacuts[n],imptimecuts[m],8e3,0,4e3,it3.GetEnergy());
                          }
                          for(int n=0;n<pidnasubcuts.size();n++){  //blob1a + blob1b (two subblobs PID_all blob1) <= there are bumps shown in the reak Mg32 886keV peak;
                            if(fbeta->fImplant.Inside(pidnasubcuts[n])){
                              FillHistogram("gamma_Na32",pidnasubcuts[n],imptimecuts[m],8e3,0,4e3,it3.GetEnergy());
                            }
                          }
                        }
                      }
                      if(m==1 &&  imptimecuts[m]->IsInside(it1.fImplantTime, it3.GetEnergy())){ // m==0 ==> Ne32
                        for(int n=0;n<pidnecuts.size();n++){ //blob1 + blob2
                          if(fbeta->fImplant.Inside(pidnecuts[n])){
                            FillHistogram("gamma_Ne30",pidnecuts[n],imptimecuts[m],8e3,0,4e3,it3.GetEnergy());
                          }
                        }
                      }
                    }
                    //if(cflag>0) continue;
                    //FillHistogram("addback", pidcuts[y], timecuts[z],8000,0,4000,it->second.AddbackE());
                    //FillHistogram("addback_implanttime", pidcuts[y], timecuts[z], 1000,0,500,it1.fImplantTime,8000,0,4000,it->second.AddbackE());
                    //cflag++;
                  }
                }             
              } 
            }
          }
        }
      }
    }
    if((x%20000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }

  }

  printf("   on entry %lu / %lu   \n",x,n);
  SaveHistograms("gamma.root");


}













