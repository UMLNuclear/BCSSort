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

//TChain *gChain = new TChain("dchan");
//TChain *gChain = new TChain("tree");
TChain *gChain = new TChain("event");
//TChain *gChain = new TChain("beta");
//TChain *gChain = new TChain("implant");
//TChain *gChain = new TChain("decay");
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
  //nentries = 1e5;
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
      if(hit.GetNumber()>=160 && hit.GetNumber()<=175){
        if(hit.GetEnergy()<600 || hit.GetNumber()==160 || hit.GetNumber()==175)
          continue;
      }
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
    //num = num.substr(0,4);
    //SaveHistograms(Form("output_%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
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
  //if(gChain->GetCurrentFile()) {
  std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());
  //================= TOFCorrection Paras Read================//
  //runnum = runnum.substr(0,4);
  //tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  //=========================================================//
  OutputManager::Get()->Set(runnum);
  //OutputManager::Get()->Set(GetRunNumber(gChain->GetCurrentFile()->GetName()));
  //std::cout << "created output tree file: " << OutputManager::Get()->GetName() << std::endl;
  std::cout << "created output tree file: " << runnum << std::endl;
  //}
  //if(tofpara.empty()) return;
  for(x=0;x<nentries;x++) {
    gChain->GetEntry(x);
    for(size_t y=0;y<event->GetNEvents();y++) {
      chan = event->GetData()[y];          
      int address      = chan->GetAddress();
      int number       = chan->GetNumber();
      double timestamp = chan->GetCoarseTime();
      double charge    = chan->GetEnergy();
      //if(number==177){
      //charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
      //}
      OutputManager::Get()->FillList(address, number, timestamp, charge);
      if(number==181){
        FillHistogram("PIN1E",10e3,0,10e3,charge);
        if(charge>100) FillHistogram("PIN1E2", 10e3,0,10e3,charge);
      }

      //FillHistogram("Summary",16000,0,16000,hit->GetEnergy(),300,0,300,hit->GetNumber());
      //FillHistogram("Summary_raw",16000,0,16000,hit->GetCharge(),300,0,300,hit->GetNumber());
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu            \r", x,nentries);
      fflush(stdout);
    } 
  }
  printf("  on entry %lu / %lu            \n", x,nentries);


  //if(gList && gROOT->GetListOfFiles()->GetSize()>0) {
  //  std::string num = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  //  num = num.substr(0,4);
  //  SaveHistograms(Form("list_output%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
  //}
  SaveHistograms(Form("listoutput%s.root",runnum.c_str()));
  OutputManager::Get()->Close();
}

//========================== Event File Sort ===========================//
// Fill implant tree and decay tree from event//
void BCSint::EventSort(){

  std::vector<TCutG*> veccut;
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/promptAll_imp0048.root");
  TCutG *_cut0 = (TCutG *)cutf2->Get("_cut0");
  TCutG *_cut1 = (TCutG *)cutf2->Get("_cut1");
  veccut.push_back(_cut0); // TOF vs timedif for 44S;
  veccut.push_back(_cut1); // Pin1E vs timedif for S;


  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent); 
  TChannel::ReadDetMapFile();

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;



  long n = gChain->GetEntries();
  //n = 1e2;
  long x = 0;

  double dt;
  for(x=0;x<n;x++) {
    gChain->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->DSSDloT()>0){ 
        dt = fevent->DSSDloT() - fevent->Pin1T();
        dt = dt/1000.;
        if(veccut[0]->IsInside(dt, fevent->I2S()) && veccut[1]->IsInside(dt, fevent->Pin1E())){
          OutputManager::Get()->FillImp(&(fevent->fHits));
        }
      }
    }else{
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        dt = fevent->HGFMax().GetTimestamp() - fevent->HGBMax().GetTimestamp();
        dt = dt/1000.;
        //if((dt>=0.04 && dt<=0.14)){ // prompt gate for decay
        if((dt>=-0.25 && dt<=-0.2) || (dt>=0.3 && dt<=0.35)){ // delay gate for decay
          OutputManager::Get()->FillDec(&(fevent->fHits));
        }
      }
    }


    if((x%200000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  OutputManager::Get()->Close();

  return;
}

