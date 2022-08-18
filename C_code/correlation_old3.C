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
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pidcut104.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  TChannel::ReadDetMapFile();
  chimp->Add("data/5us_tofcor/sample/imp_good/implant104*");
  chdec->Add("data/5us_tofcor/sample/dec_prompt/decay104*");

  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);

  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 1e3;
  //double minT = -1;

  int f = 11;
  int b = 19;
  int cflag = 0;  
  int c = 0;
  int decsize = 0;  

  pixel piximp;

  ofstream ofile;
  ofile.open("correlation.txt");
  

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
      cutname = "Na";
    }else if(cut1[1]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "Ne";
    }else{
      cutname = "unknown";
    }
    //ofile<<"isotope = "<< cutname<<std::endl;
    //ofile<<"Imp(" << ximp-1 <<") at " <<forimp.DSSDloT()<<std::endl;; //ximp is not entry of "forimp";
    ofile<<"Implant["<< forimp.GetPixel().first<<","<<forimp.GetPixel().second<<"] at "<< forimp.DSSDloT()<<" ";

    //printf("isotope = %s\n", cutname.c_str());
    //printf("Imp(%lu) at %lu\n", ximp-1, (size_t)forimp.DSSDloT()); //ximp is not entry of "forimp";
// get decay: didn't check decay come earlier than "forimp" or not, it should be fine    
  decsize=0;
  while(xdec<ndec){
      chdec->GetEntry(xdec++);
      if(fdec->DSSDhiT() > fimp->DSSDloT()) {xdec--; break;}
      if(fdec->DSSDhiT() < forimp.DSSDloT()) continue;
      if(fdec->GetPixel().first!=f || fdec->GetPixel().second!=b) continue;
      decsize++;
      double sumE = 0;
      for(auto &it:fdec->fGe) {if(it.GetEnergy()>10 && it.GetEnergy()<4000)sumE += it.GetEnergy();} 
      //printf("\t\t\tDec(%lu) at %lu(dt=%f) with %zu gamma(%fkeV)\n",xdec-1, (size_t)fdec->DSSDhiT(), (fdec->DSSDhiT()-forimp.DSSDloT())/1e6,fdec->fGe.size(),sumE);
      //ofile<<"\t\t\tDec("<< xdec-1 <<") at "<<fdec->DSSDhiT() <<"(dt="<< (fdec->DSSDhiT()-forimp.DSSDloT())/1e6<<") with "<< fdec->fGe.size()<<" gamma("<<sumE<<"keV)"<<std::endl;
    }
    fimp->Copy(forimp);
    //printf("\n");
    ofile<<decsize<<" decay"<<std::endl;
    if(c>5) break; 
  }
  ofile.close();
  return;
}

