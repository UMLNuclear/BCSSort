#include<map>
#include<cstdio>

#include<TChain.h>
#include<TFile.h>
#include<TTree.h>
#include<TCutG.h>
#include<TKey.h>
#include<TList.h>

#include<TChannel.h>
#include<Implant.h>
#include<TOFCorrection.h>
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

std::map<pixel,Implant*> fImpMap;
std::map<pixel,std::vector<Decay> *> fDecMap;
std::map<pixel, double> fImpTime;

//==========Initial Map=========//
void Initial(){
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new Implant;
      fDecMap[pix] = new std::vector<Decay>;
      fImpTime[pix] = -1;
    }
  }
}

//==========correlation=========//
void correlation(){
  int runnum = 1040;
  std::vector<double> tofpara = TOFCorrection::Get()->ReadFile(runnum, "/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
  if(tofpara.empty()){
    printf("run%i.root file doesn't have tof parameters\n", runnum);
    return;
  } 
  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pidcut104.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/prompt_Imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); // TOF vs timedif for 44S;
  veccut.push_back(pin1E_dt); // Pin1E vs timedif for S;

  TChannel::ReadDetMapFile();
  chimp->Add(Form("data/5us_tofcor/implant/imp_good/implant%i*",runnum));
  chdec->Add(Form("data/5us_tofcor/decay/dec_prompt/decay%i*",runnum));
  //chdec->Add(Form("data/5us_tofcor/sample/decay/dec_delay/decay%i*",runnum));
  
  TFile *newf = new TFile(Form("beta_good_prompt%i.root",runnum),"recreate");
  //TFile *newf = new TFile(Form("beta_good_delay%i.root",runnum),"recreate");

  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);

  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 5e3;
  //double minT = -1;
  double dt;

  TTree *tree = new TTree("beta","beta tree");
  Beta fbeta;
  tree->Branch("Beta",&fbeta);

  Initial();

  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    fimp->fI2S = fimp->fI2S + tofpara[0];
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
          fbeta.Set(*current, *vdec);
          dt = fbeta.fImplant.DSSDloT() - fbeta.fImplant.fPIN1T;
          dt = dt/1000.; //us
          if(veccut[0]->IsInside(dt, fbeta.fImplant.fI2S) && veccut[1]->IsInside(dt,fbeta.fImplant.fPIN1E)){
            tree->Fill();
          }
          current->Clear();
          vdec->clear();
          break;
        }
        for(int i=pixdec.first-1;i<=pixdec.first+1;i++){
          if(i<0 || i>39) continue;
          for(int j=pixdec.second-1;j<=pixdec.second+1;j++){
            if(j<0 || j>39) continue;
            pixel temp = std::make_pair(i,j);
            if(fImpTime[temp]>0 && fImpTime[temp]<fdec->DSSDhiT()){
              fDecMap[temp]->push_back(*fdec);
            }
          }
        }
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
      fbeta.Set(*current, *vdec);
      dt = fbeta.fImplant.DSSDloT() - fbeta.fImplant.fPIN1T;
      dt = dt/1000.; //us
      if(veccut[0]->IsInside(dt, fbeta.fImplant.fI2S) && veccut[1]->IsInside(dt,fbeta.fImplant.fPIN1E)){
        tree->Fill();
      }
      fbeta.Clear();
      current->Clear();
      vdec->clear();
    }
  } 

  printf("  on entry %lu / %lu  \n", ximp, nimp);

  tree->Write();
  newf->Close();
  return;
}

