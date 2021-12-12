#include<map>
#include<cstdio>

#include<TChain.h>
#include<TFile.h>
#include<TCutG.h>
#include<TKey.h>
#include<TList.h>

#include<TChannel.h>
#include<Implant.h>
#include<globals.h>
#include<util.h>


// global variables;
TChain *chimp = new TChain("implant");
TChain *chdec = new TChain("decay");

Implant *fimp = new Implant;
Decay *fdec = new Decay;

long nimp, ndec;
long ximp = 0;
long xdec = 0;



//==========correlation=========//
void correlation_old3(){

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("root_file/cuts/pidcut_48.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  TChannel::ReadDetMapFile();
  chimp->Add("run0048/10us/no160175_Esssd600/imp_good/implant0048-0*");
  //chdec->Add("run0048/10us/no160175_Esssd600/dec_prompt/decay0048-0*");
  chdec->Add("run0048/10us/no160175_Esssd600/dec_delay/decay0048-0*");

  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);

  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 1e3;
  //double minT = -1;

  int f = 0;
  int b = 4;
  int cflag = 0;  
  int c = 0;

  pixel piximp;

//get first implant at "piximp";
  while(piximp.first!=f || piximp.second!=b){
    chimp->GetEntry(ximp++);
    piximp = fimp->GetPixel();
  }
  c++;
  Implant forimp(*fimp); 
// get the second implant;
  while(ximp<nimp){
    chimp->GetEntry(ximp++);
    if(fimp->GetPixel().first!=f || fimp->GetPixel().second!=b) continue;
    c++;
    
    string cutname;
    if(cut1[0]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "42P";
    }else if(cut1[1]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "43S";
    }else if(cut1[2]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "44S";
    }else if(cut1[3]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "45Cl";
    }else{
      cutname = "unknown";
    }
    printf("isotope = %s\n", cutname.c_str());
    printf("Imp(%lu) at %lu\n", ximp-1, (size_t)forimp.DSSDloT()); //ximp is not entry of "forimp";
// get decay: didn't check decay come earlier than "forimp" or not, it should be fine    
  while(xdec<ndec){
      chdec->GetEntry(xdec++);
      if(fdec->DSSDhiT() > fimp->DSSDloT()) {xdec--; break;}
      if(fdec->DSSDhiT() < forimp.DSSDloT()) continue;
      if(fdec->GetPixel().first!=f || fdec->GetPixel().second!=b) continue;
      double sumE = 0;
      for(auto &it:fdec->fGe) {if(it.GetEnergy()>10 && it.GetEnergy()<4000)sumE += it.GetEnergy();} 
      printf("\t\t\tDec(%lu) at %lu with %zu gamma(%fkeV)\n",xdec-1, (size_t)fdec->DSSDhiT(), fdec->fGe.size(),sumE);
    }
    fimp->Copy(forimp);
    printf("\n");
    if(c>50) break; 
  }
  return;
}

