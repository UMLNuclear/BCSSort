#include<map>
#include<cstdio>

#include<TChain.h>
#include<TFile.h>
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

std::map<pixel, Implant*> fImpMap;
std::map<pixel, std::vector<Decay> *> fDecMap1;
std::map<pixel, std::vector<Decay> *> fDecMap3;
std::map<pixel, double> fTimeMap; // record implant pixel is implanted or not; former implant time;

//=======Initialize maps========//
void Initial(){
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new Implant;
      fDecMap1[pix] = new std::vector<Decay>;
      fDecMap3[pix] = new std::vector<Decay>;
      fTimeMap[pix] = -1;
    }
  }
}

//============= Match implant and decays===========//
void Match(double minT, double maxT){
  for(xdec=0;xdec<ndec;xdec++){
    double dt = DBL_MAX;
    pixel best = std::make_pair(-1,-1);
    chdec->GetEntry(xdec);
    if(fdec->DSSDhiT()<minT) continue;
    if(fdec->DSSDhiT()>maxT) return;
    pixel pixdec = fdec->GetPixel(); // front: 0~39 * back: 0~39;
    for(int row=pixdec.first-1;row<=pixdec.first+1;row++){
      for(int col=pixdec.second-1;col<=pixdec.second+1;col++){
        if((row<0 || row>39) || (col<0 || col>39)) continue;
        pixel temp = std::make_pair(row, col);
        if(fTimeMap[temp]>0){ // that pixel is implanted
          if(fdec->DSSDhiT()>fTimeMap[temp]){ // decay must come later than implant
            if(temp.first==pixdec.first && temp.second==pixdec.second) fDecMap1[temp]->push_back(*fdec);
            if((fdec->DSSDhiT() - fTimeMap[temp])<dt){ // get the time closest implant;
              dt = fdec->DSSDhiT() - fTimeMap[temp];
              best = temp;
            }       
          }
        }
      }
    }
    if(best.first<0 || best.second<0){ // decay doesn't have correlated implant;
    }else{
      fDecMap3[best]->push_back(*fdec);
    }
  }
}

//==========correlation=========//
void correlation_old2(){

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("root_file/cuts/pidcut_48.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  TChannel::ReadDetMapFile();
  chimp->Add("run0048/10us/no160175_Esssd600/imp_good/implant0048-0*");
  chdec->Add("run0048/10us/no160175_Esssd600/dec_prompt/decay0048-0*");

  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);

  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();
  //nimp = 1e3;
  //double minT = -1;

  int f = 12;
  int b = 15;
  int cflag = 0;  
  int c = 0;
  Initial();


  pixel piximp;
  while(piximp.first!=f || piximp.second!=b){
    chimp->GetEntry(ximp++);
    piximp = fimp->GetPixel();
  }
  c++;
  Implant forimp(*fimp); 
  while(ximp<nimp){
    chimp->GetEntry(ximp++);
    if(fimp->GetPixel().first!=f || fimp->GetPixel().second!=b) continue;
    c++;
    
    string cutname;
    if(cut1[0]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "42P";
    }else if(cut1[1]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "43S";
    }else if(cut1[2]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "44S";
    }else if(cut1[3]->IsInside(forimp.fI2S, forimp.fPIN1E)) {
      cutname = "45Cl";
    }else{
      cutname = "unknown";
    }
    printf("isotope = %s\n", cutname.c_str());
    printf("Imp(%lu) at %lu\n", ximp-2, (size_t)forimp.DSSDloT());

    while(xdec<ndec){
      chdec->GetEntry(xdec++);
      if(fdec->DSSDhiT() > fimp->DSSDloT()) {xdec--; break;}
      if(fdec->GetPixel().first!=f || fdec->GetPixel().second!=b) continue;
      printf("\t\t\tDec(%lu) at %lu\n",xdec-1, (size_t)fdec->DSSDhiT());
    }
    fimp->Copy(forimp);
    printf("\n");
    if(c>5) break; 
  }
  return;
}

