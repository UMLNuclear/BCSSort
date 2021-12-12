#include "Implant.h"
#include"TChannel.h"
#include<TH2D.h>
#include<TH2.h>
#include<TDirectory.h>
#include<TCutG.h>
#include<TFile.h>

#include<cstdio>
#include<map>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <algorithm>

#include <util.h>

Implant::Implant() {
    Clear();
}
Implant::~Implant() {}

Implant::Implant(std::vector<DetHit> *hits){
    Set(hits);
}


////////////////////////////////////////////////////////
/*HPGeHit::HPGeHit(DetHit hit) {
  GetNumber()    = hit.GetNumber();
  fAddress   = hit.fAddress;
  fCharge    = hit.fCharge; 
  GetTimestamp() = hit.GetTimestamp();
  }
  HPGeHit::~HPGeHit() {}
  */



bool Implant::Stopped() const {
  bool stopped = true;
  for(auto it : fSSSD) {
    if(it.GetCharge()>100)
      return false;
  }
  return stopped;
}





////// Set (buill event from <vector>DetHit) ////////
void Implant::Set(std::vector<DetHit> *hits){
    //DetHit strip;
    for (int y=0; y<hits->size(); y++){// get each element from event
        DetHit  hitt = hits->at(y);
        switch (hitt.GetNumber()){
            case 40 ... 79:
                fDSSDFront.push_back(DetHit(hitt));
                break;

            case 120 ... 159:
                fDSSDBack.push_back(DetHit(hitt));
                break;

            case 160 ... 175:
                fSSSD.push_back(DetHit(hitt));
                break;

            case 176:
                fI2N = hitt.GetCharge();
                fI2NT = hitt.GetTimestamp();
                break;
            
            case 177:
                fI2S = hitt.GetCharge();
                fI2ST = hitt.GetTimestamp();
                break;

            case 180:
                fI2S_I2N = hitt.GetCharge();
                fI2S_I2N_T = hitt.GetTimestamp();
                break;

            case 181:
                fPIN1E = hitt.GetCharge();
                fPIN1T = hitt.GetTimestamp();
                break;

            case 182:
                fPIN2E = hitt.GetCharge();
                fPIN2T = hitt.GetTimestamp();
                break;

            default:
                break;
        } 
    }
}




//// Clear (clear all variables) //////////
void Implant::Clear(){
    fI2S = 0;
    fI2ST = 0;
    fI2N = 0;
    fI2NT = 0;
    fPIN1E = 0;
    fPIN1T = 0;
    fPIN2E = 0;
    fPIN2T = 0;
    fI2S_I2N = 0;
    fI2S_I2N_T = 0;
    fDSSDFront.clear();
    fDSSDBack.clear();
    fSSSD.clear(); 
}


void Implant::SimplePrint() const{
    printf("(%02i: %02i) Eng[%.1f,%.1f]  @ %f\n", fDSSDFront[0].GetNumber()-40, fDSSDBack[0].GetNumber()-120, fDSSDFront[0].GetEnergy(), fDSSDBack[0].GetEnergy(),fDSSDFront[0].GetTimestamp() );
}


////////// Print /////////////
void Implant::Print() const{
    if(IsGood()) {
     std::cout << "GOOD IMPLANT:\t"; 
    } else {
     std::cout << "BAD IMPLANT:\t"; 
    }
    printf("FS = %i BS = %i pin1 @ %f \n",FrontSize(),BackSize(),GetTimestamp());
    if(IsGood()) {
      printf("\t(%02i: %02i) @ %.1f \t FHE = %f\t BHE = %f \n",
              GetPixel().first,
              GetPixel().second, 
              //fDSSDFront[0].GetNumber(), 
              //fDSSDBack[0].GetNumber()-80, 
              fDSSDFront[0].GetTimestamp(), 
              fDSSDFront[0].GetCharge(), 
              fDSSDBack[0].GetCharge()); 
    } else {
    }

}



void Implant::PrintFun(TCutG *mycut) const {
    if(!mycut) printf("no TCutG\n");
    else{
        if(mycut->IsInside(fPIN1E,fPIN2E)){
            printf("implant in %s \n", mycut->GetName());
        }   
        else printf("implant is not in %s\n", mycut->GetName());
    }
    printf("------------------------\n");
    printf("Implant @ %f\n",fPIN1T);
    printf("\tpin1 energy:\t%f\n",fPIN1E);
    printf("\tpin2 energy:\t%f\n",fPIN2E);
    printf("\t%lu front strips fired\n",fDSSDFront.size());  
    printf("\t%lu back strips fired\n",fDSSDBack.size());  
    int x=0;
    if(fDSSDFront.size()>fDSSDBack.size()){x=fDSSDFront.size();}
    else x=fDSSDBack.size();
    for(int i=0;i<x;i++){
        printf("%i:\t",i);
        if(i<fDSSDFront.size()){
            printf("%02i:  %.1f\t",fDSSDFront[i].GetNumber(), fDSSDFront[i].GetEnergy());
        }else{ 
            printf("\t\t\t");
        }
        if(i<fDSSDBack.size()){
            printf("%02i:  %.1f\n",fDSSDBack[i].GetNumber(), fDSSDBack[i].GetEnergy());
        }else{ 
            printf("\n");
        }
    }
    printf("------------------------\n");
}


//////// Draw (return one TH2D)////////////
TH2 *Implant::Draw(){

    TH2 *fl= (TH2D*)gDirectory->FindObject("BL_FL_flE");
    if(!fl){
        fl = new TH2D("BL_FL_flE","BL_FL_flE", 40,0,39, 40,0,39);
    }
    fl->Reset();
    fl->SetTitle(Form("DSSD@ %lu", fPIN1E));
    //TH2D *fl = new TH2D("BL_FL_flE", "BL_FL_flE", 40,0,39, 40,0,39);
    //TH2D *bl = new TH2D("BL_FL_blE", "BL_FL_blE", 40,0,39, 40,0,39);   
    for (int i=0;i<BackSize();i++){
        for (int j=0; j<FrontSize();j++){
            fl->Fill(fDSSDBack[i].GetNumber()-120, fDSSDFront[j].GetNumber()-40, fDSSDFront[j].GetCharge());
            //bl->Fill(stripb[i].GetNumber()-120, stripf[j].GetNumber()-40);
        }
    }

    return fl;   
}

bool Implant::IsGood() const { // for 44S, PID + prompt gate*2
    bool good = false;
    static std::vector<TCutG*> veccut;
    if(veccut.size()==0) {
      TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/pidcut_48.root");
      TCutG * _cut2 = (TCutG *)cutf1->Get("_cut2"); 
      veccut.push_back(_cut2);// PID for 44S; 
      TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/promptCut_imp0048.root");
      TCutG * _cut0 = (TCutG *)cutf2->Get("_cut0");
      TCutG * _cut1 = (TCutG *)cutf2->Get("_cut1");
      veccut.push_back(_cut0); // TOF vs timedif for 44S;
      veccut.push_back(_cut1); // Pin1E vs timedif for S;
    }
    if(veccut[0]->IsInside(fI2S,fPIN1E)){
      if(DSSDloT()>0){
        double dt = DSSDloT() - fPIN1T;
        dt = dt/1000.;
        if((veccut[1]->IsInside(dt, fI2S)) && veccut[2]->IsInside(dt, fPIN1E)){
          good = true;
        }
      }  
    }

    return good;
}

double Implant::DSSDloT() const {
  double t = -1;
  if(FrontSize()>0){
    t = LGFMax().GetTimestamp();
  }else if(BackSize()>0){
    t = LGBMax().GetTimestamp();
  }
  return t;
}

DetHit Implant::LGFMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:fDSSDFront){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

DetHit Implant::LGBMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:fDSSDBack){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

std::pair<int, int> Implant::GetPixel() const{
    //printf(__PRETTY_FUNCTION__ );
    //printf("\n");
    int indexf = -1;
    int indexb = -1;
    if(FrontSize()>0) indexf = -LGFMax().GetNumber()+79;
    if(BackSize()>0) indexb  = -LGBMax().GetNumber()+159;
    return std::make_pair(indexf,indexb);
}


void Implant::Copy(Implant &lhs) const{
    lhs.fI2S   = this->fI2S  ;
    lhs.fI2ST  = this->fI2ST ;
    lhs.fI2N   = this->fI2N  ;
    lhs.fI2NT  = this->fI2NT ;
    lhs.fPIN1E = this->fPIN1E;
    lhs.fPIN1T = this->fPIN1T;
    lhs.fPIN2E = this->fPIN2E;
    lhs.fPIN2T = this->fPIN2T;
    lhs.fI2S_I2N = this->fI2S_I2N;
    lhs.fI2S_I2N_T = this->fI2S_I2N_T;
    lhs.fDSSDFront = this->fDSSDFront;
    lhs.fDSSDBack  = this->fDSSDBack ;
    lhs.fSSSD      = this->fSSSD     ; 
}









// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
Decay::Decay(){fImplantTime = -1;}
Decay::~Decay(){}



////// Set (buill event from <vector>DetHit) ////////
void Decay::Set(std::vector<DetHit> *hits){
    DetHit strip;
    for (int y=0; y<hits->size(); y++){// get each element from event
        DetHit  hitt = hits->at(y);
        switch (hitt.GetNumber()){
            case 0 ... 39:
                fDSSDFront.push_back(DetHit(hitt)); 
                break;

            case 80 ... 119:
                fDSSDBack.push_back(DetHit(hitt)); 
                break;

            case 160 ... 175:
                fSSSD.push_back(DetHit(hitt)); 
                break;

            case 192 ... 207:
                fLaBr.push_back(DetHit(hitt));
                break;

            case 208 ... 271:
                fGe.push_back(DetHit(hitt));
                break;

            default:
                break;
        } 
    }
}



//// Clear (clear all variables) //////////
void Decay::Clear(){
    fDSSDFront.clear();
    fDSSDBack.clear();
    fLaBr.clear();
    fGe.clear();
    fSSSD.clear(); 
}


bool Decay::IsPrompt() const{
    bool good = false;    
    if(FrontSize()>0 && BackSize()>0){
      double dt = HGFMax().GetTimestamp() - HGBMax().GetTimestamp();
      dt = dt/1000.;
      if(dt>=0.04 && dt<=0.14){
        good = true;
      }
    }

    return good;
}

bool Decay::IsDelay() const{
    bool good = false;    
    if(FrontSize()>0 && BackSize()>0){
      double dt = HGFMax().GetTimestamp() - HGBMax().GetTimestamp();
      dt = dt/1000.;
      if(fabs(dt+0.225)<=0.025 || fabs(dt-0.325)<=0.025){
        good = true;
      }
    }

    return good;
}

double Decay::DSSDhiT() const {
  double t = -1;
  if(FrontSize()>0){
    t = HGFMax().GetTimestamp();
  }else if(BackSize()>0){
    t = HGBMax().GetTimestamp();
  }
  return t;
}

DetHit Decay::HGFMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:fDSSDFront){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

DetHit Decay::HGBMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:fDSSDBack){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

std::pair<int,int> Decay::GetPixel() const{
    int indexf = -1;
    int indexb = -1;
    if(FrontSize()>0) indexf = HGFMax().GetNumber();
    if(BackSize()>0) indexb  = HGBMax().GetNumber()-80;
    return std::make_pair(indexf,indexb); 
}


void Decay::SimplePrint(double dt) const{
    printf("(%02i: %02i) @ %.1f \t FHE = %f\t BHE = %f \t dt(ms) = %.1f\n", fDSSDFront[0].GetNumber(), fDSSDBack[0].GetNumber()-80, fDSSDFront[0].GetTimestamp(), fDSSDFront[0].GetEnergy(), fDSSDBack[0].GetEnergy(), dt*1e-6 );
}


// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //
// ============================================================= //

//Clover::Clover(){}
//Clover::~Clover(){}

std::vector<double> Clover::AddbackSum(TCutG *gate) const{
    double dt = 0;
    double eref=0;
    int n = 0;
    double sum = 0;
    std::vector<double> energy;
    for (size_t x=0;x<Size();x++){
        energy.push_back(fXtal[x].GetEnergy());
    }
    for(size_t x=0;x<Size();x++){
        for(size_t y=x+1;y<Size();y++){
            if(fXtal[x].GetEnergy()>fXtal[y].GetEnergy()){
                dt = fXtal[x].GetTimestamp() - fXtal[y].GetTimestamp();
                eref = fXtal[y].GetEnergy();
            }else{
                dt = fXtal[y].GetTimestamp() - fXtal[x].GetTimestamp();
                eref = fXtal[x].GetEnergy();
            }
            if(gate->IsInside(dt,eref)){
                n+=1;
                sum += fXtal[x].GetEnergy();
                sum += fXtal[y].GetEnergy();
                energy[x] = 0;
                energy[y] = 0;
            }
        }
    }
    n = n/3+1;
    sum = sum/n;
    energy.erase(std::remove(energy.begin(), energy.end(),0), energy.end());
    energy.push_back(sum);
    return energy; 
}




/*int Clover::ReadXtalFile(const char* filename, Option_t *opt){

  std::string infilename = filename;
  std::fstream infile;
  std::string line;
  infile.open(infilename.c_str());
  int i=0;
  while(getline(infile,line)){
  int det;
  double c0;
  double c1;
  double c2;
  double c3;

  std::stringstream ss(line);
  ss >> det;
  ss >> c0;
  ss >> c1;
  ss >> c2;
  ss >> c3;

  xtalmat[det].push_back(c0);
  xtalmat[det].push_back(c1);
  xtalmat[det].push_back(c2);
  xtalmat[det].push_back(c3);
  i +=1;
  }

  return i;
  }*/



// --------------------------------------------------------------//
// --------------------------------------------------------------//
// --------------------------------------------------------------//
// --------------------------------------------------------------//
// --------------------------------------------------------------//
// --------------------------------------------------------------//


void Beta::Clear(){
    fImplant.Clear();
    for(int x=0;x<DecaySize();x++){
        fDecay[x].Clear();
    }
}

void Beta::SimplePrint() const{
    printf("======================\n");
    printf("implant:\n");
    fImplant.SimplePrint();
    for(size_t x=0;x<DecaySize();x++){
        printf("\t");
        fDecay[x].SimplePrint(fDecay[x].fDSSDFront[0].GetTimestamp()-fImplant.fDSSDFront[0].GetTimestamp());
    }
}






