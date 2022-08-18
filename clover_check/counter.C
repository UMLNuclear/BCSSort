

#include<util.h>
#include<TChain.h>
#include<map>
#include<cstdio>







void counter(){

  TChannel::ReadDetMapFile();
  
  TChain *chlist = new TChain("tree");
  chlist->Add("clover_check/data_5us/list1040-00.root");
  int fnumber;
  double fcharge;
  chlist->SetBranchAddress("number", &fnumber);
  chlist->SetBranchAddress("charge", &fcharge);
  long n_list = chlist->GetEntries();

  TChain *cheve = new TChain("event");
  cheve->Add("clover_check/data_5us/event1040-00.root");
  BCSEvent *fevent = new BCSEvent;
  cheve->SetBranchAddress("BCSEvent", &fevent);
  long n_event = cheve->GetEntries();

  long x = 0;
  long c1 = 0;
  long c2 = 0;
  long c3 = 0;
  long c4 = 0;
  long c5 = 0;
  long c6 = 0;
  long c7 = 0;
  long c8 = 0;
  long c9 = 0;
  long c10 = 0;
  long c11 = 0;
  long c12 = 0;
  long c13 = 0;
  long c14 = 0;
  long c15 = 0;
  long c16 = 0;
  long c17 = 0;
  long c18 = 0;
  long c19 = 0;
  long c20 = 0;

  std::map<int,std::vector<double>> clmap;
  
  for(x=0;x<n_list;x++){
    chlist->GetEntry(x);
    if(fnumber>207 && fnumber<272){
      clmap[fnumber].push_back(fcharge);
    }
  }
  
  printf("channel \t counts\n");
  std::map<int,std::vector<double>>::iterator it;
  for(it=clmap.begin();it!=clmap.end();it++){
    printf("%i \t %lu\n", it->first, it->second.size());
  }


}
