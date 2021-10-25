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

  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  long x=0;
  long n = gChain->GetEntries();
  cout<<n<<endl;

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    Implant fimp = fbeta->fImplant;
    FillHistogram("PIN1E", 10e3,0,10e3,fimp.fPIN1E);
    if(fimp.fPIN1E>100) FillHistogram("PIN1E2", 10e3,0,10e3,fimp.fPIN1E);
    if(fimp.Stopped()){
      FillHistogram("PIN1E_nosssd", 10e3,0,10e3,fimp.fPIN1E);
      if(fimp.fPIN1E>100) FillHistogram("PIN1E_nosssd2", 10e3,0,10e3,fimp.fPIN1E);
      if(fimp.IsGood()){
        FillHistogram("PIN1E_isgood_nosssd", 10e3,0,10e3,fimp.fPIN1E);
        if(fimp.fPIN1E>100) FillHistogram("PIN1E_isgood_nosssd2", 10e3,0,10e3,fimp.fPIN1E);
      }
    }

    if((x%20000)==0) {
      printf("on entry %lu / %lu   \r",x,n);
      fflush(stdout);
    }

  }

  printf("   on entry %lu / %lu   \n",x,n);
  SaveHistograms("betaoutput1031.root");



}
