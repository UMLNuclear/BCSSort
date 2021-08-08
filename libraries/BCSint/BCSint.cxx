#include <pwd.h>
#include <string>

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

#include <globals.h>

//TChain *gChain = new TChain("dchan");
TChain *gChain = new TChain("event");
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
  //nentries = 1e6;
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
      //std::cout << "===  Start ===" << std::endl;
      FillHistogram("Summary",16000,0,16000,hit.GetEnergy(),300,0,300,hit.GetNumber());
      FillHistogram("Summary_raw",16000,0,16000,hit.GetCharge(),300,0,300,hit.GetNumber());
      //hit.print();
      //std::cout << "dt = " << hit.GetTimestamp() << " - " << starttime << " = " <<  hit.GetTimestamp()-starttime <<std::endl << std::endl;
      //if(  hit.GetTimestamp()-starttime < 0) exit(0);
      if((hit.GetTimestamp()-starttime)>buildtime){
        Correlator::Get()->AddEvent(vec_hit);
        starttime = hit.GetTimestamp();
        vec_hit->clear();
        //std::cout << " End Correlator->AddEvent " << std::endl << std::endl;
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
    SaveHistograms(Form("output%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
  }

  OutputManager::Get()->Close();
  //    CloseCoorMap();
}

void BCSint::TOFfluctuation(){
  //outfile.append(GetRunNumber(gChain->GetCurrentFile()->GetName()));
  
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  std::string num = runnum.substr(0,4);  

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long x = 0;
  long n = gChain->GetEntries();
  double first_time = -1;
  double runtime;
  double tof;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
        runtime = (fevent->Pin1T()-first_time)/1e9;
        tof = fevent->I2S();
        FillHistogram(Form("tof%s",num.c_str()), 3800,0,3800,runtime, 1600,0,32000,tof);
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("  on entry %lu / %lu  \n",x,n);
  //SaveHistograms("tof.root","update");
  SaveHistograms(Form("tof%s.root", num.c_str()));

}

void BCSint::CorrectTOF(){
  
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  std::string num = runnum.substr(0,4);  
  
  TFile *corf = TFile::Open(Form("tof%s.root",num.c_str()));
  TH2D *hist = (TH2D *)corf->Get(Form("tof%s",num.c_str()));
  if(hist==NULL){
    printf("histogram is empty\n");
    return;
  }
  TProfile *pfx = hist->ProfileX();
  TSpline3 *sp3 = new TSpline3(pfx);
  TF1 px("px","pol1");
  pfx->Fit(&px);
  double offset = px.Eval(0);

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long x = 0;
  long n = gChain->GetEntries();
  double first_time = -1;
  double runtime;
  double ctof,tof;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(first_time<0) first_time = fevent->Pin1T();
    if(fevent->I2S()>0){
        runtime = (fevent->Pin1T()-first_time)/1e9;
        tof = fevent->I2S();
        ctof = tof + (offset-sp3->Eval(runtime));
        FillHistogram(Form("ctof%s", num.c_str()), 3600,0,3600,runtime, 3000,0,30000,ctof);
    }
    if((x%50000)==0){
      printf("  on entry %lu / %lu \r",x,n);
      fflush(stdout);
    }
  }

  printf("  on entry %lu / %lu  \n",x,n);
  //SaveHistograms("ctof.root","update");
  SaveHistograms(Form("ctof%s.root",num.c_str()));

}














