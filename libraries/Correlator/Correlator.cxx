#include <string>
#include <cstdio>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cfloat>

#include<TFile.h>
#include<TChain.h>
#include<TTree.h>
#include<TCutG.h>
#include<TKey.h>
#include<TList.h>

#include<TChannel.h>
#include <BCSint.h>
#include <Correlator.h>
#include <DetHit.h>
#include <Implant.h>
#include <util.h>
#include <OutputManager.h>

#include <globals.h>

Correlator *Correlator::fCorrelator=0;

Correlator *Correlator::Get() { 
  if(!fCorrelator) {
    fCorrelator = new Correlator;
  }
  return fCorrelator;
}

Correlator::Correlator() {
  fImplant = new Implant;
  fDecay   = new Decay;
  InitMaps();
}

Correlator::~Correlator() { }


void Correlator::AddEvent(std::vector<DetHit> *event) {
  OutputManager::Get()->FillEvent(event);
  bool pin1     = false;
  bool lofront  = false;
  bool loback   = false;
  bool hifront  = false;
  bool hiback   = false;
  bool sssdlo   = false;
  bool clarray  = false;
  double pin1E = -1;
  double i2s = -1;
  int cpin1=0;
  for(auto &it : *event) {
    FillHistogram("sum_event_raw", 20000,0,40000,it.GetCharge(), 300,0,300,it.GetNumber());
    FillHistogram("sum_event_cal", 16000,0,8000,it.GetEnergy(), 300,0,300,it.GetNumber());
    if(it.GetNumber()==181) {
      cpin1++;
      pin1E = it.GetCharge();
    }else if(it.GetNumber()==177){
      i2s = it.GetCharge();
    } 
  }
  if(pin1E>0){
    FillHistogram("counts_pin1",100,0,100,cpin1);
    if(i2s>0)
      FillHistogram("pid",2e3,0,2e4,i2s,8e3,0,8e3,pin1E);
  }

  fImplant->Clear();
  fDecay->Clear();
  return;
}

void Correlator::AddImplant(std::vector<DetHit> *event) { 
  fImplant->Set(event);
  pixel pix = fImplant->GetPixel(); // Get low gain pixel with max energy; [0,39]*[0,39]
  if((pix.first>=0 && pix.first<40) && (pix.second>=0 && pix.second<40)){ // check this pixel is "goodPos" (valid)
    if(fImplantMap[pix]->IsGood()){ //this pixel is implanted. imp at pix must fire both LGF and LGB;
      double timeDiffIon = fImplant->fPIN1T - fImplantMap[pix]->fPIN1T; // time difference between two imp.
      if(timeDiffIon>0){ // the current must come later tha the existing imp; // back-to-back imp.
        if(timeDiffIon>minImplantDT) FlushPixel(pix,true);
        fBlockMap[pix] = timeDiffIon; 
      }
    }else{// if this pixle is empty
      fImplant->Copy(*fImplantMap[pix]);
    }
  }
  return;
}

void Correlator::AddDecay(std::vector<DetHit> *event) { 
  fDecay->Set(event);
  pixel pix = Match(fDecay->GetPixel(),fDecay->GetTimestamp());
  if(pix.first<0 || pix.second<0) { return; } //failed to find match.
  fDecay->SetDecayTime( (fDecay->DSSDhiT() -  fImplantMap[pix]->GetTimestamp())   );
  fDecayMap[pix]->push_back(*fDecay);
  return;
}

void Correlator::AddCloverOnly(std::vector<DetHit> *event){
  fDecay->Set(event);
  double cloverTime = event->at(0).GetTimestamp(); // first hit time as current_time;
  double isomerTimeCut = 40e3; //unit:ns = 40us;
  for(auto &it:fImplantMap){
    if(it.second->IsGood()){
      double isomerTDiff = cloverTime - it.second->GetTimestamp();
      if(isomerTDiff>0 && isomerTDiff<isomerTimeCut){
        fDecay->SetDecayTime(isomerTDiff); 
        fDecayMap[it.first]->push_back(*fDecay);  
      } 
    }
  } 
  return; 
}

pixel Correlator::Match(pixel dpix,double T) {  //given a decay, returns mathcing implant.
  double dt= DBL_MAX;
  pixel best = std::make_pair(-1,-1);
  for(int row=dpix.first-1;row<=dpix.first+1;row++) {
    for(int col=dpix.second-1;col<=dpix.second+1;col++) {
      if((row<0||row>39) || (col<0||col>39)) continue;
      pixel temp = std::make_pair(row,col);
      if(fImplantMap[temp]->IsGood()){ // decay must came later than imp;
        double ionDecayTDiff = fImplantMap[temp]->fPIN1T - T;
        if(ionDecayTDiff>0){ 
          if(ionDecayTDiff <dt) {
            dt = ionDecayTDiff;
            best = temp;
          }
        }
      }
    }
  }
  // check this pixel again
  if(best.first<0 || best.first>39 || best.second<0 || best.second>39) best =  std::make_pair(-1,-1);
  else{
    if(dt<0) best = std::make_pair(-1,-1); // if decay comes earlier than imp;
    else if(dt > tDiffCorrCut) best = std::make_pair(-1,-1); // if decay delays too long;
    else if(fBlockMap[best]>0 && fBlockMap[best]<minImplantDT) best = std::make_pair(-1,-1); // if the imp is back-to-back
  }
  return best;
}

void Correlator::InitMaps() {
  for(int x=0;x<40;x++){
    for(int y=0;y<40;y++){
      pixel pix = std::make_pair(x,y);
      fImplantMap[pix] = new Implant;
      fDecayMap[pix]   = new std::vector<Decay>;    
      fBlockMap[pix]   = -1;
    }
  }    
}



void Correlator::CleanImplantMap(double dT) {
  for(auto &it : fImplantMap) {
    if(!it.second->IsGood()) continue;
    double tdiff = fabs(dT-it.second->fPIN1T);
    if(tdiff>EXPIRETIME) { //remove(clear) Implant from the map
      Implant *current = it.second;
      FlushPixel(it.first,true);
    }
  }
  for(auto &it : fBlockMap) {
    if(it.second<0) continue; 
    double tdiff = fabs(dT-it.second);
    if(tdiff>EXPIRETIME) it.second = -1;
  }
}

void Correlator::FlushPixel(pixel pix, bool histogram) {
  Implant            *current = fImplantMap[pix];
  std::vector<Decay> *vdec    = fDecayMap[pix];  

  if(histogram) {
    //=======================================//
    OutputManager::Get()->Fill(current,vdec);
    //=======================================//
    //FillHistogram("pid_nogate", 2e3,0,2e4,current->fI2S, 5e3,0,20e3,current->fPIN1E);
    //for(auto &it1:*vdec){
    //  for(auto &it2:it1.fGe){
    //    FillHistogram("single_imptime", 3000,0,1500,it1.fImplantTime,8000,0,4000,it2.GetEnergy());
    //  }
    //}
    for(size_t m=0;m<vdec->size();m++){
      //if(vdec->at(m).FrontSize()==0) continue;
      for(auto &it:vdec->at(m).fGe){
        FillHistogram("single_imptime", 3000,0,1500,vdec->at(m).fDecayTime,8000,0,4000,it.GetEnergy());
      }
    }

    //=======================================//
    //=======================================//
  }
  current->Clear();
  vdec->clear();
}

void Correlator::FlushAll(bool histogram) {
  for(int x=0;x<40;x++) { 
    for(int y=0;y<40;y++) { 
      pixel pix = std::make_pair(x,y);
      if(fImplantMap[pix]->IsGood()) FlushPixel(pix,histogram);
    }
  }
  return;
}



