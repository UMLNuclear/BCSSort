#include <string>
#include <cstdio>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cfloat>

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
  for(auto &it1 : *event){
    FillHistogram("chan_sum", 8e3,0,16e3,it1.GetCharge(), 300,0,300,it1.GetNumber());
  }
  bool pin1=false;
  for(auto &it : *event) {
    if(it.GetNumber()==181) {  
      pin1=true; continue;// change break to continue
    }
  }
  if(pin1) { 
    AddImplant(event);
  } else { 
    AddDecay(event);
  }
  fImplant->Clear();
  fDecay->Clear();
  return;
}

void Correlator::AddImplant(std::vector<DetHit> *event) { 
  fImplant->Set(event);
  CleanImplantMap(fImplant->fPIN1T);
  if(fImplant->IsGood()) {
    pixel pix = fImplant->GetPixel();
    //CleanImplantMap(fImplant->fPIN1T);
    if(fImplantMap[pix]->IsGood() || (fBlockMap[pix]>0) ) { 
      fBlockMap[pix] = fImplant->fPIN1T; 
      FlushPixel(pix,false);
      return;
    }
    fImplant->Copy(*fImplantMap[pix]);
  }
  return;
}

void Correlator::AddDecay(std::vector<DetHit> *event) { 
  fDecay->Set(event);
  if(fDecay->IsGood()) {
    pixel pix = Match(fDecay->GetPixel(),fDecay->GetTimestamp());
    if(pix.first<0 || pix.second<0) { return; } //failed to find match.
    fDecay->SetImplantTime( (fDecay->GetTimestamp() -  fImplantMap[pix]->GetTimestamp())/1e6   );
    fDecayMap[pix]->push_back(*fDecay);


  }

  return;
}

pixel Correlator::Match(pixel dpix,double T) {  //given a decay, returns mathcing implant.
  double dt= DBL_MAX;
  pixel best = std::make_pair(-1,-1);
  int count=0;
  for(int row=dpix.first-1;row<=dpix.first+1;row++) {
    for(int col=dpix.second-1;col<=dpix.second+1;col++) {
      if((row<0||row>39) || (col<0||col>39)) continue;
      pixel temp = std::make_pair(row,col);
      if(fImplantMap[temp]->IsGood() && fImplantMap[temp]->Stopped()) {
        count++;
        if(fabs(fImplantMap[temp]->GetTimestamp()-T) <dt) {
          dt = fImplantMap[temp]->GetTimestamp()-T;
          best = temp;
        }
      }
    }
  }
  FillHistogram("count",10,0,10,count);
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
    FillHistogram("pid_cur", 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    if(current->Stopped()){
      FillHistogram("pidnosssd_cur", 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    }


    TIter iter(gCuts);
    while(TCutG *cut = (TCutG*)iter.Next()) {
      if(!cut->IsInside(current->fI2S,current->fPIN1E)) continue;    
      FillHistogram("pid",cut, 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
      if(current->Stopped()) 
        FillHistogram("pidnosssd",cut, 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
      for(auto &it1 : *vdec) {

         FillHistogram("summary",cut, 300,0,300,it1.fDSSDFront[0].GetNumber(), 8000,0,4000,it1.fDSSDFront[0].GetEnergy());
         FillHistogram("summary_raw",cut, 300,0,300,it1.fDSSDFront[0].GetNumber(), 8000,0,64000,it1.fDSSDFront[0].GetCharge());
         FillHistogram("summary",cut, 300,0,300,it1.fDSSDBack[0].GetNumber(), 8000,0,4000,it1.fDSSDBack[0].GetEnergy());
         FillHistogram("summary_raw",cut, 300,0,300,it1.fDSSDBack[0].GetNumber(), 8000,0,64000,it1.fDSSDBack[0].GetCharge());


        for(auto &it2 : it1.fGe) { 
          double dt = it1.GetTimestamp() - it2.GetTimestamp();

          FillHistogram("timedif_decay_gamma",cut,200,-1000,1000,dt);
          FillHistogram("timedif2_decay_gamma",cut,200,-1000,1000,dt,4000,0,4000,it2.GetEnergy());

          FillHistogram("summary",cut, 300,0,300,it2.GetNumber(), 8000,0,4000,it2.GetEnergy());
          FillHistogram("summary_raw",cut, 300,0,300,it2.GetNumber(), 8000,0,64000,it2.GetCharge());
 
          //double deltatime = fabs(it1.GetTimestamp() - it2.GetTimestamp());

          //if(deltatime>40 && deltatime<300){ // Fill gamma spectrum
          if(dt>-200 && dt<300){ // Fill gamma spectrum
            FillHistogram("summary_prompt",cut, 300,0,300,it2.GetNumber(), 8000,0,4000,it2.GetEnergy());
            FillHistogram("summary_prompt_raw",cut, 300,0,300,it2.GetNumber(), 8000,0,64000,it2.GetCharge());
          }
        }
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
