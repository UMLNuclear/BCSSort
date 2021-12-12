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

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("root_file/cuts/pidcut_48.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  TChannel::ReadDetMapFile();
  chimp->Add("run0048/10us/no160175_Esssd600/imp_good/implant0048-0*");
  //chdec->Add("run0048/10us/no160175_Esssd600/dec_prompt/decay0048-0*");
  chdec->Add("run0048/10us/no160175_Esssd600/dec_delay/decay0048-0*");

  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);

  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 5e3;
  //double minT = -1;

  //TFile *newf = new TFile("beta_prompt0048.root","recreate");
  TFile *newf = new TFile("beta_delay0048.root","recreate");
  TTree *tree = new TTree("beta","beta tree");
  Beta fbeta;
  tree->Branch("Beta",&fbeta);

  Initial();

  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    pixel piximp = fimp->GetPixel();
    if(piximp.first<0 ||  piximp.second<0) continue;
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
          tree->Fill();
          //if(piximp.first==12 && piximp.second==15) printf("Imp[%lu] at %lu \n", ximp, (size_t)fImpMap[piximp]->DSSDloT());
          fbeta.Clear();
          current->Clear();
          vdec->clear();
          break;
        }
        if(fdec->DSSDhiT() < fImpTime[pixdec]) continue; // dceay must come later than implant;
        fDecMap[pixdec]->push_back(*fdec);
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
      tree->Fill();
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

