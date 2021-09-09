#include <pwd.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSpline.h>

#include <DDASEvent.h>
#include <BCSint.h>
#include <BCSintFunctions.h>
#include <DetHit.h>
#include <TChannel.h>
#include <Correlator.h>
#include <util.h>
#include <OutputManager.h>
#include <TOFCorrection.h>

#include <globals.h>

TChain *gChain = new TChain("hit_tree");
//TChain *gChain = new TChain("dchan");
//TChain *gChain = new TChain("event");
//TChain *gChain = new TChain("beta");
TList  *gCuts  = new TList();


BCSint *BCSint::fInstance=0;

BCSint *BCSint::Get() {
  if(fInstance==0) 
    fInstance = new BCSint();
  return fInstance;
}

BCSint::BCSint(const char *app,int *argc,char **argv,void *options,int numOptions,bool noLogo) 
  : TRint(app,argc,argv,options,numOptions,noLogo) { 

    fFileCount=0;

    SetPrompt();

  }   

BCSint::~BCSint() { }


const char *BCSint::SetPrompt(const char *newPrompt) {

  return TRint::SetPrompt(newPrompt);
}


void BCSint::Terminate(int status) {


  //Be polite when you leave.
  printf("\nbye,bye\t%s\n",getpwuid(getuid())->pw_name);


  TRint::Terminate(status);

}

//void BCSint::OpenRootFile(std::string fname) {
//  std::string rfile = Form("_file%i",fFileCount++);
//  //gROOT->ProcessLine(Form("TFile *%s = new TFile(\"%s\");",rfile.c_str(),fname.c_str()));  
//  gROOT->ProcessLine(Form("TFile *%s = TFile::Open(\"%s\");",rfile.c_str(),fname.c_str()));  
//  std::cout << "file " << fname << " opened as " << rfile << std::endl;
//  //TFile *f = (TFile*)gROOT->FindObjectAny(rfile.c_str());
//  TFile *f = (TFile*)gROOT->GetListOfFiles()->Last();
//  //std::cout << f->GetName() << std::endl;
//  if(f->FindObjectAny("dchan")) {
//    gChain->Add(fname.c_str());
//  }
//}


void BCSint::DoSort() {

  //if(!gChain || gChain->GetEntries()==0) return;

  //    InitCoorMap();
  LoadCuts();
  DDASEvent *event  = new DDASEvent;
  ddaschannel *chan = new ddaschannel;

  gChain->SetBranchAddress("ddasevent",&event);

  long nentries = gChain->GetEntries();
  nentries = 1e3;
  long x = 0;

  if(gChain->GetCurrentFile()) {
    OutputManager::Get()->Set(GetRunNumber(gChain->GetCurrentFile()->GetName()));
    std::cout << "created output tree file: " << OutputManager::Get()->GetName() << std::endl;
  }

  double starttime = 0;
  double buildtime = BUILDTIME;
  std::vector<DetHit> *vec_hit = new std::vector<DetHit>;

  TChannel::ReadDetMapFile();
  for(x=0;x<nentries;x++) {
    gChain->GetEntry(x);
    for(size_t y=0;y<event->GetNEvents();y++) {
      chan = event->GetData()[y];          
      DetHit hit(chan);
      FillHistogram("Summary",16000,0,16000,hit.GetEnergy(),300,0,300,hit.GetNumber());
      FillHistogram("Summary_raw",16000,0,16000,hit.GetCharge(),300,0,300,hit.GetNumber());
      if((hit.GetTimestamp()-starttime)>buildtime){
        Correlator::Get()->AddEvent(vec_hit);
        starttime = hit.GetTimestamp();
        vec_hit->clear();
      }
      vec_hit->push_back(hit);
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu            \r", x,nentries);
      fflush(stdout);
    } 
  }
  Correlator::Get()->FlushAll();
  printf("  on entry %lu / %lu            \n", x,nentries);


  if(gList && gROOT->GetListOfFiles()->GetSize()>0) {
    std::string num = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
    num = num.substr(0,4);
    SaveHistograms(Form("output%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
  }

  OutputManager::Get()->Close();
  //    CloseCoorMap();
}


//=======================Fill List.file After TOFCorrection================//
void BCSint::ListSort(){ 
  LoadCuts();
  DDASEvent *event  = new DDASEvent;
  ddaschannel *chan = new ddaschannel;

  gChain->SetBranchAddress("ddasevent",&event);

  long nentries = gChain->GetEntries();
  //nentries = 1e3;
  long x = 0;

  TChannel::ReadDetMapFile();
  std::vector<double> tofpara;
  if(gChain->GetCurrentFile()) {
    std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());
//================= TOFCorrection Paras Read================//
    runnum = runnum.substr(0,4);
    tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
//=========================================================//
    OutputManager::Get()->Set(runnum);
    //OutputManager::Get()->Set(GetRunNumber(gChain->GetCurrentFile()->GetName()));
    //std::cout << "created output tree file: " << OutputManager::Get()->GetName() << std::endl;
    std::cout << "created output tree file: " << runnum << std::endl;
  }
  //printf("%f : %f :%f\n", tofpara[0], tofpara[1], tofpara[2]);
  if(tofpara.empty()) return;
  for(x=0;x<nentries;x++) {
    gChain->GetEntry(x);
    for(size_t y=0;y<event->GetNEvents();y++) {
      chan = event->GetData()[y];          
      DetHit *hit = new DetHit(chan);
      if(hit->GetNumber()==177){
        //printf("old tof = %f\n",hit->GetCharge());
        double i2s = hit->GetCharge();
        i2s = tofpara[0]+tofpara[1]*i2s+tofpara[2]*i2s*i2s;
        hit->SetCharge(i2s);
        //printf("new tof = %f  \n\n",hit->GetCharge());
      }
      //if(hit->GetNumber()==180) printf("I2S_I2N = %f\n\n", hit->GetCharge());
      OutputManager::Get()->FillList(hit);
         
      FillHistogram("Summary",16000,0,16000,hit->GetEnergy(),300,0,300,hit->GetNumber());
      FillHistogram("Summary_raw",16000,0,16000,hit->GetCharge(),300,0,300,hit->GetNumber());
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu            \r", x,nentries);
      fflush(stdout);
    } 
  }
  printf("  on entry %lu / %lu            \n", x,nentries);


  if(gList && gROOT->GetListOfFiles()->GetSize()>0) {
    std::string num = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
    num = num.substr(0,4);
    SaveHistograms(Form("list_output%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
  }
  OutputManager::Get()->Close();
}













