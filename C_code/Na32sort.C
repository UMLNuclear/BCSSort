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
#include<TOFCorrection.h>

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
  vector<TCutG*> cuts7;
  //TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/blobcuts.root"); //cuts from pid
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newsubblobNa.root"); //cuts from pid, after TOF fluctuation and TOF&I2 correction;
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
  TFile *mycut6 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/PID_Ne30_cuts.root"); // cuts from pid only Na30 
  TIter keys6 (mycut6->GetListOfKeys());
  while(TKey *key6 = (TKey*)keys6.Next()) {
    pidnecuts.push_back((TCutG*)key6->ReadObj());
  }
  TFile *mycut7 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newblobcuts.root"); // cuts from newpid after mom&TOF correction
  TIter keys7 (mycut7->GetListOfKeys());
  while(TKey *key7 = (TKey*)keys7.Next()) {
    cuts7.push_back((TCutG*)key7->ReadObj());
  }

  //n=1e5;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    int pidflag = 0;
    int pidflag1 = 0;
    if(fbeta->fImplant.Stopped()){
      for(auto &it1: fbeta->fDecay){
        int i=-1;
        for(auto &it2: it1.fGe){
          int j=-1; //count for ggmat
          i++;
          double dt = it1.GetTimestamp() - it2.GetTimestamp();
          if(timecuts[0]->IsInside(dt, it2.GetEnergy()) && fbeta->fImplant.Inside(pidnacuts[0])){ //LHE + blob1(only 32Na)
            //if(it1.fImplantTime<100){
              double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
              dtime = fabs(dtime)/1e6; //unit:ms
              FillHistogram("decaytime",1000,0,500,dtime);
              FillHistogram("single_implanttime", 1000,0,500,it1.fImplantTime,8000,0,4000,it2.GetEnergy());
              for(auto &it3 : it1.fGe){
                j++;
                if(i<j){
                  FillHistogram("ggmatrix", 4000,0,4000,it2.GetEnergy(),4000,0,4000,it3.GetEnergy());
                  FillHistogram("ggmatrix", 4000,0,4000,it3.GetEnergy(),4000,0,4000,it2.GetEnergy());
                }
              }
              if(pidflag1>0) continue;
              FillHistogram("pid_1",imptimecuts[0], 2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
              pidflag1++;
            //}
          }
          //FillHistogram("single_implanttime",timecuts[0],1000,0,500,it1.fImplantTime,8000,0,4000,it2.GetEnergy());
          if(timecuts[0]->IsInside(dt, it2.GetEnergy()) && fbeta->fImplant.Inside(pidcuts[0])){ //LHE + blob1
            if(imptimecuts[0]->IsInside(it1.fImplantTime,it2.GetEnergy())){ //Mg32-886keV        
              if(pidflag>0) continue;
              FillHistogram("pid",imptimecuts[0], 2e3,0,2e4,fbeta->fImplant.fI2S, 8e3,0,8e3,fbeta->fImplant.fPIN1E);
              pidflag++;
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
    SaveHistograms("Na32new.root");


    }













