#ifndef __CINT__
#include <map>

#include<cstdio>
#include<string>
#include<iostream>

#include <DDASEvent.h>
#include<DetHit.h>
#include<TChannel.h>
#include<util.h>

#include<TChain.h>
#include<TFile.h>
#endif



void ddashcheck(int argc=1, const char** argv=0){



  TChannel::ReadDetMapFile();

  TChain *chain = new TChain("dchan");
  if(argc<1){
    for(int i=1;i<argc;i++){
      chain->Add(argv[i]);
    }
  }else{
    chain->Add("/home/zhu/packages/BCSSort/CORAnalysis/data/run1045-00.root");
  }
  DDASEvent *event  = new DDASEvent;
  ddaschannel *chan = new ddaschannel;

  chain->SetBranchAddress("ddasevent",&event);
  long nentries = chain->GetEntries();
  //if(nentries==0) return 0;
  long x = 0;
  double time = 0;

  std::multimap<double,DetHit> mm;

  while(x<nentries){
    chain->GetEntry(x++);
    for(size_t y=0;y<event->GetNEvents();y++) {
      chan = event->GetData()[y];
      DetHit hit(chan);
      double time = chan->GetCoarseTime();
      mm.insert(std::make_pair(time,hit));
      //if((hit.GetTimestamp()-time)<0){
      //  printf("entry = %lu @ %zu \t\t current = %f(former = %f) \n", x,y,hit.GetTimestamp(),time);
      //  printf("in ddashcnanel:   \t\t current = %f \n\n", chan->GetCoarseTime());
      //}
      //time = hit.GetTimestamp();
  
    }

    //if((x%5000)==0){
    //  printf(" on entry = %lu / %lu\r",x,nentries);
    //  fflush(stdout);
    //}


  }

  for(auto it: mm) {
    if((it.second.GetTimestamp()-time)<0){
      //printf("entry = %lu @ %zu \t\t current = %f(former = %f) \n", x,y,hit.GetTimestamp(),time);
      //printf("in ddashcnanel:   \t\t current = %f \n\n", chan->GetCoarseTime());
      printf("BAD\n");
    }

    time = it.second.GetTimestamp();
  }



  printf(" on entry = %lu / %lu\n",x,nentries);


}


#ifndef __CINT__

int main(int argc, const char** argv) {
  ddashcheck(argc, argv);
  return 0;
}


#endif



