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
#include <Implant.h>
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
  std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());
  runnum = runnum.substr(0,4);
  std::vector<double> tofpara;
  //tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  //if(tofpara.empty()){
  //  printf("run%i.root file doesn't have tof parameters\n", runnum);
  //  return;
  //} 

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
      FillHistogram("sum_das_raw",20000,0,40000,chan->GetEnergy(),300,0,300,chan->GetNumber());
      int address      = chan->GetAddress();
      int number       = chan->GetNumber();
      double timestamp = chan->GetCoarseTime();
      double charge    = chan->Energy();
      //if(number==177) {
      //  FillHistogram("tof_list_before",2e3,0,2e4, charge);
      //  charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
      //  FillHistogram("tof_list_after",2e3,0,2e4, charge);
      //}
      OutputManager::Get()->FillList(address, number, timestamp, charge);
      DetHit hit(chan);
      if(hit.GetNumber()==177){
        FillHistogram("tof_event_before",2e3,0,2e4, hit.GetCharge());
        double c = hit.GetCharge();
        c = tofpara[0]+tofpara[1]*c+tofpara[2]*c*c;
        hit.SetCharge(c);
        FillHistogram("tof_event_after",2e3,0,2e4, hit.GetCharge());
      }
      if((hit.GetTimestamp()-starttime)>buildtime){
        Correlator::Get()->AddEvent(vec_hit);
        starttime = hit.GetTimestamp();
        vec_hit->clear();
      }
      vec_hit->push_back(hit);
    }
    if((x%200000)==0){
      printf("  on entry %lu / %lu            \r", x,nentries);
      fflush(stdout);
    } 
  }
  Correlator::Get()->FlushAll();
  printf("  on entry %lu / %lu            \n", x,nentries);


  if(gList && gROOT->GetListOfFiles()->GetSize()>0) {
    std::string num = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
    num = num.substr(0,4);
    //SaveHistograms(Form("tof_output_%s.root",num.c_str())); // saves any histograms to output.root, will over write EVER time.
  }

  OutputManager::Get()->Close();
}


//========================== Event File Sort ===========================//
// Fill implant tree and decay tree from event//
void BCSint::EventSort(){

  std::vector<TCutG*> veccut;
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/prompt_Imp104.root");
  TCutG *tof_dt = (TCutG *)cutf2->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf2->Get("pin1E_dt");
  veccut.push_back(tof_dt); 
  veccut.push_back(pin1E_dt); 


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
        //OutputManager::Get()->FillImp(&(fevent->fHits));
        dt = fevent->DSSDloT() - fevent->Pin1T();
        dt = dt/1000.;
        //FillHistogram("Pin1E_dt", 1000,-5,5,dt,10e2,0,10e3,fevent->Pin1E());
        //FillHistogram("TOF_dt", 1000,-5,5,dt,20e2,0,20e3,fevent->I2S());
        if(veccut[0]->IsInside(dt, fevent->I2S()) && veccut[1]->IsInside(dt, fevent->Pin1E())){
          OutputManager::Get()->FillImp(&(fevent->fHits));
        }
      }
    }else{
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        //OutputManager::Get()->FillDec(&(fevent->fHits));
        dt = fevent->HGFMax().GetTimestamp() - fevent->HGBMax().GetTimestamp();
        dt = dt/1000.;
        FillHistogram("dec_dt_hiFandB", 1000,-5,5,dt);
        //if((dt>=-0.22 && dt<=-0.11) || (dt>=0.29 && dt<=0.4)){ // delay gate for decay
        if(dt>=0.03 && dt<=0.14){ // prompt gate for decay
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
  //SaveHistograms(Form("event5op_dt_%s.root",runnum.c_str()));

  return;
}

//================ Correlation from imp_good and dec_prompt=============//
//===gChain = "event"===//
void BCSint::Correlation(){
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;

  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/prompt_Imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); // TOF vs timedif for 44S;
  veccut.push_back(pin1E_dt); // Pin1E vs timedif for S;

  std::map<pixel, Implant*> fImpMap;
  std::map<pixel, std::vector<Decay> *> fDecMap;
  std::map<pixel, double> fImpTime;

  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new Implant;
      fDecMap[pix] = new std::vector<Decay>;
      fImpTime[pix] = -1;
    }
  }

  TChannel::ReadDetMapFile();
  TChain *chdec = new TChain("decay");
  TChain *chimp = new TChain("implant");
  chimp->Add(Form("data/5us_tofcor/implant/implant%s*",runnum.c_str()));
  std::cout<< "Add implant file: "<< runnum << std::endl;
  chdec->Add(Form("data/5us_tofcor/decay/dec_prompt/decay%s*",runnum.c_str()));
  std::cout<< "Add decay file: "<< runnum << std::endl;
  Implant *fimp = new Implant;
  Decay *fdec = new Decay;
  chimp->SetBranchAddress("Implant",&fimp);
  chdec->SetBranchAddress("Decay",&fdec);

  long nimp,ndec;
  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 1e3; 
  long ximp = 0;
  long xdec = 0;
  double dt;

  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    pixel piximp = fimp->GetPixel();
    if(piximp.first<0 ||  piximp.second<0) continue;
    if(piximp.first>39 ||  piximp.second>39) continue;
    if(fImpTime[piximp]>0){ //that pixel is implanted;
      while(xdec<ndec){
        chdec->GetEntry(xdec++);
        pixel pixdec = fdec->GetPixel();
        if(pixdec.first<0 || pixdec.second<0) continue;
        if(fdec->DSSDhiT() > fimp->DSSDloT()) {
          xdec--;
          Implant *current = fImpMap[piximp];
          std::vector<Decay> *vdec = fDecMap[piximp];
          dt = current->DSSDloT() - current->fPIN1T;
          dt = dt/1000.; //us
          if(veccut[0]->IsInside(dt, current->fI2S) && veccut[1]->IsInside(dt,current->fPIN1E)){
            OutputManager::Get()->Fill(current,vdec);
          }
          current->Clear();
          vdec->clear();
          break;
        }
        double tshort = DBL_MAX;
        pixel best = std::make_pair(-1,-1);
        for(int i=pixdec.first-1;i<=pixdec.first+1;i++){
          if(i<0 || i>39) continue;
          for(int j=pixdec.second-1;j<=pixdec.second+1;j++){
            if(j<0 || j>39) continue;
            pixel temp = std::make_pair(i,j);
            if(fImpTime[temp]>0 && fImpTime[temp]<fdec->DSSDhiT()){
              if(fabs(fImpTime[temp] - fdec->DSSDhiT())<tshort){
                tshort = fabs(fImpTime[temp] - fdec->DSSDhiT());
                best = temp;
              }
            }
          }
        }
        if(best.first<0 || best.first>39) continue;
        if(best.second<0 || best.second>39) continue;
        dt = fdec->DSSDhiT() - fImpMap[pixdec]->fPIN1T;
        dt = dt/1000000.;//unit:ms 
        fdec->SetDecayTime(dt);
        fDecMap[best]->push_back(*fdec);
        
        /*for(int i=pixdec.first-1;i<=pixdec.first+1;i++){
          if(i<0 || i>39) continue;
          for(int j=pixdec.second-1;j<=pixdec.second+1;j++){
            if(j<0 || j>39) continue;
            pixel temp = std::make_pair(i,j);
            if(fImpTime[temp]>0 && fImpTime[temp]<fdec->DSSDhiT()){
              dt = fdec->DSSDhiT() - fImpMap[pixdec]->fPIN1T;
              dt = dt/1000000.;//unit:ms 
              fdec->SetDecayTime(dt);
              fDecMap[temp]->push_back(*fdec);
            }
          }
        }*/



      }
    }
    fimp->Copy(*fImpMap[piximp]);
    fImpTime[piximp] = fimp->DSSDloT();

    if((ximp%5000)==0){
      printf("  on entry %lu / %lu    \r", ximp, nimp);
      fflush(stdout);
    }
  } 

  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      if(fImpTime[pix]<0) continue;
      Implant *current = fImpMap[pix];
      std::vector<Decay> *vdec = fDecMap[pix];
      dt = current->DSSDloT() - current->fPIN1T;
      dt = dt/1000.; //us
      if(veccut[0]->IsInside(dt, current->fI2S) && veccut[1]->IsInside(dt,current->fPIN1E)){
        OutputManager::Get()->Fill(current,vdec);
      }
      current->Clear();
      vdec->clear();
    }
  }

  printf("  on entry %lu / %lu  \n", ximp, nimp);
  OutputManager::Get()->Close();

  return;
}

//============================== from beta tree to beta tree by beta_tof correction ==========================//
void BCSint::BetaSort(){

  TChannel::ReadDetMapFile();
  
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  TChain *beta = new TChain("beta");
  beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt/correlation1/beta%s*.root",runnum.c_str()));
  Beta *fbeta = new Beta;
  beta->SetBranchAddress("Beta",&fbeta);
  std::vector<double> tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum), 
                                                                "/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
  if(tofpara.empty()){
    printf("run%i.root file doesn't have tof parameters\n", runnum);
    return;
  } 
  OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;
  
  long x = 0;
  long n = beta->GetEntries();
  
  for(x=0;x<n;x++){
    beta->GetEntry(x);
    Implant *fimp = &(fbeta->fImplant);
    FillHistogram("pid_before", 4000,0,20000,fimp->fI2S, 4000,0,8000,fimp->fPIN1E);
    fimp->fI2S = fimp->fI2S + tofpara[0];
    FillHistogram("pid_after", 4000,0,20000,fimp->fI2S, 4000,0,8000,fimp->fPIN1E);
    OutputManager::Get()->Fill(fimp, &(fbeta->fDecay));

    if((x%5000)==0){
      printf("  on entry %lu / %lu    \r", x, n);
      fflush(stdout);
    }  
  }
  printf("  on entry %lu / %lu  \n", x, n);
  OutputManager::Get()->Close();
  SaveHistograms(Form("beta_op%s.root",runnum.c_str()));  


}

















