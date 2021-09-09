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




void Na32sort (){

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
  TFile *mycut2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/decaytimecuts.root"); // cuts from decaytime and HPGe.Energy; timecuts[0] is "LHE" cut
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
    int pidflag = 0;
    int pidflag1 = 0;
    int pidflag2 = 0;
    if(fbeta->fImplant.Stopped()){
      //========= blobcuts from pid============//
      //for(long y=0;y<pidcuts.size();y++){
      //if(fbeta->fImplant.Inside(pidcuts[y])){
      for(auto &it1: fbeta->fDecay){
        for(auto &it2: it1.fGe){
          double dt = it1.GetTimestamp() - it2.GetTimestamp();
          if(imptimecuts[0]->IsInside(it1.fImplantTime,it2.GetEnergy())){ //Mg32-886keV
            if(!pidflag1){
              FillHistogram("pid_1",imptimecuts[0], 2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
              pidflag1++;
            }
            if(imptimecuts[0]->IsInside(it1.fImplantTime,it2.GetEnergy()) && timecuts[0]->IsInside(dt, it2.GetEnergy())){ //Mg32-886keV + LHE
              if(!pidflag2){
                FillHistogram("pid_2",imptimecuts[0], 2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
                pidflag2++;
              }
            }
            //FillHistogram("single_implanttime",timecuts[0],1000,0,500,it1.fImplantTime,8000,0,4000,it2.GetEnergy());
            if(pidflag>0) continue;
            if(timecuts[0]->IsInside(dt, it2.GetEnergy()) && fbeta->fImplant.Inside(pidcuts[0])){ //LHE + blob1
              if(imptimecuts[0]->IsInside(it1.fImplantTime,it2.GetEnergy())){ //Mg32-886keV
                FillHistogram("pid",imptimecuts[0], 2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
                pidflag++;
              }
            }
          } 
        }
        //}
    }

    }
    if((x%20000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }

    }

    printf("   on entry %lu / %lu   \n",x,n);
    SaveHistograms("Na32.root");


  }













