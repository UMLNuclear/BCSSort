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

long nimp, ndec, ximp, xdec;

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
void correlation(){
  
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
  //nimp = 1e4;
  double minT = -1;

  int f = 25;
  int b = 25;
  int cflag = 0;  

  Initial();

  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    pixel piximp = fimp->GetPixel(); // front:0~39 * back:0~39;

    if(fTimeMap[piximp]>0){ // check that pixel is empty or not;
      Match(minT, fimp->DSSDloT());
      minT = fimp->DSSDloT();
      std::vector<Decay> *vdec1 = fDecMap1[piximp];
      std::vector<Decay> *vdec3 = fDecMap3[piximp]; 

      if(piximp.first==f && piximp.second==b){
        cflag ++;
        double dtimp = fimp->DSSDloT() - fImpMap[piximp]->DSSDloT();
        dtimp = dtimp/1000000.;
        printf("================= Implant ========================\n");
        if(cut1[0]->IsInside(fImpMap[piximp]->fI2S, fImpMap[piximp]->fPIN1E)) {
          printf("Entry = Implant(42P) at [%i,%i] %luns exits for %fms\n",piximp.first, piximp.second, fImpMap[piximp]->DSSDloT(),dtimp);
        }else if(cut1[1]->IsInside(fImpMap[piximp]->fI2S, fImpMap[piximp]->fPIN1E)) {
          printf("Implant(43S) at [%i,%i] exits for %fms\n",piximp.first, piximp.second, dtimp);
        }else if(cut1[2]->IsInside(fImpMap[piximp]->fI2S, fImpMap[piximp]->fPIN1E)) {
          printf("Implant(44S) at [%i,%i] exits for %fms\n",piximp.first, piximp.second, dtimp);
        }else if(cut1[3]->IsInside(fImpMap[piximp]->fI2S, fImpMap[piximp]->fPIN1E)) {
          printf("Implant(45Cl) at [%i,%i] exits for %fms\n",piximp.first, piximp.second, dtimp);
        }else{
          printf("Implant(nothing) at [%i,%i] exits for %fms\n",piximp.first, piximp.second, dtimp);
        }
        printf("================= Decay at 1*1 ====================\n");
        printf("Size = %zu\n", vdec1->size());
        for(int c=0;c<vdec1->size();c++){
          double dtdec = vdec1->at(c).DSSDhiT() - fImpMap[piximp]->DSSDloT();
          dtdec = dtdec/1000000.;
          pixel pixdec = vdec1->at(c).GetPixel();
          printf("Decay[%i] at [%i,%i] %luns delay %fms\n",c, pixdec.first, pixdec.second, vdec1->at(c).DSSDhiT(),dtdec);
        }
        printf("================= Decay at 3*3 ====================\n");
        printf("Size = %zu\n", vdec3->size());
        for(int c=0;c<vdec3->size();c++){
          double dtdec = vdec3->at(c).DSSDhiT() - fImpMap[piximp]->DSSDloT();
          dtdec = dtdec/1000000.;
          pixel pixdec = vdec3->at(c).GetPixel();
          printf("Decay[%i] at [%i,%i] %luns delay %fms\n",c, pixdec.first, pixdec.second, vdec3->at(c).DSSDhiT(),dtdec);
        }
        printf("\n");
      }
      if(cflag >= 1) return;
      vdec1->clear();
      vdec3->clear();
    }
    fimp->Copy(*fImpMap[piximp]);
    fTimeMap[piximp] = fimp->DSSDloT();

  } 
  return;
}

