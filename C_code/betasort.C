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






void betasort(){

  TChannel::ReadDetMapFile();
  int runnum = 1042;
  
  TChain *beta2 = new TChain("beta");
  beta2->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_tofcor/correlation2bestT/beta%i*.root",runnum));
  Beta *fbeta2 = new Beta;
  beta2->SetBranchAddress("Beta",&fbeta2);
  long n2 = beta2->GetEntries();

  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)mycut1->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");
  


  long x=0;
  double dE, dnum, dt, sumE;
  int c0 = 0;
  int c1 = 0;
  int c2 = 0;
  int c3 = 0;


  for(x=1000;x<1100;x++){
    beta2->Clear();
    beta2->GetEntry(x);
    Implant fimp2;
    fimp2.Clear();
    fimp2 = fbeta2->fImplant;
    //pixel piximp2 = fimp2.GetPixel();
    pixel piximp2 = fbeta2->fImplant.GetPixel();
    for(int m=0;m<fbeta2->DecaySize();m++){
      bool flag2 = false;
      for(auto &it:fbeta2->fDecay[m].fGe){
        if(it.GetEnergy()>10 && it.GetEnergy()<4000){
          flag2 = true;
          break;
        }
      }
      if(!flag2) continue;
      pixel pixdec2 = fbeta2->fDecay[m].GetPixel();
      c2++;
      double decaytime2 = fbeta2->fDecay[m].fDecayTime;
      if(fabs(pixdec2.first-piximp2.first)>2 || fabs(pixdec2.second-piximp2.second)>2){
        c3++;
        printf("size = 2\nrunnum=%i\nentry=%lu\n", runnum, x);
        printf("Implant[%i,%i]\tDecay[%i][%i,%i] @ %fms\n", fimp2.GetPixel().first, fimp2.GetPixel().second, m, fbeta2->fDecay[m].GetPixel().first, fbeta2->fDecay[m].GetPixel().second, fbeta2->fDecay[m].fDecayTime);
        printf("Implant[%i,%i]\tDecay[%i][%i,%i] @ %fms\n\n", piximp2.first, piximp2.second, m, pixdec2.first, pixdec2.second, decaytime2);
      }
    }

    //if((x%5000)==0) {
    //  printf("on entry %lu / %lu   \r",x,n2);
    //  fflush(stdout);
    //}

  }

  printf("   on entry %lu / %lu   \n",x,n2);
  std::cout<<"decay = "   << c2 <<std::endl;
  std::cout<<"Unsatisfy correlation = " << c3 <<std::endl;
  //SaveHistograms(Form("check_beta_tree%i.root",runnum));

}
