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
  bool pin1     = false;
  bool lofront  = false;
  bool loback   = false;
  bool hifront  = false;
  bool hiback   = false;
  bool sssdlo   = false;
  bool clarray  = false;
  double pin1E = -1;
  for(auto &it : *event) {
    if(it.GetNumber()==181) {  
      pin1=true; 
      pin1E = it.GetEnergy();
      continue;// change break to continue
    }else if(it.GetNumber()>=0 && it.GetNumber()<40){ // dssd_hi_front
      hifront = true; continue;
    }else if(it.GetNumber()>39 && it.GetNumber()<80){ // dssd_lo_front
      lofront = true; continue;
    }else if(it.GetNumber()>79 && it.GetNumber()<120){ // dssd_hi_back
      hiback = true; continue;
    }else if(it.GetNumber()>119 && it.GetNumber()<160){ // dssd_lo_back
      loback = true; continue;
    }else if(it.GetNumber()>271 && it.GetNumber()<288){ //sssd_lo
      sssdlo = true; continue;
    }else if(it.GetNumber()>207 && it.GetNumber()<272){ // clovers
      clarray = true; continue;
    }
  }
  if(hifront || lofront || hiback || loback){ //hasDSSD 
    if(pin1 && (lofront && loback) && (pin1E>4000 && pin1E<30000)){ //isImplant
      AddImplant(event);
    }else if(!pin1 && (hifront && hiback) && !sssdlo){ //isDecay
      AddDecay(event);
    }
  }else{
    if(!pin1 && clarray){
      AddCloverOnly(event);
    }
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
  fDecay->SetImplantTime( (fDecay->GetTimestamp() -  fImplantMap[pix]->GetTimestamp())/1e6   );
  fDecayMap[pix]->push_back(*fDecay);
  return;
}

void Correlator::AddCloverOnly(std::vector<DetHit> *event){
  fDecay->Set(event);
  double cloverTime = event[0].GetTimestamp(); // first hit time as current_time;
  double isomerTimeCut = 40e3; //unit:ns = 40us;
  for(auto &it:fImplantMap){
    if(it.second->IsGood()){
      double isomerTDiff = cloverTime - it.second->GetTimestamp();
      if(isomerTDiff>0 && isomerTDiff<isomerTimeCut){
        fDecay->SetImplantTime(isomerTDiff/1e6); 
        fDecayMap[it.first]->push(*fDecay);  
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
      if(ionDecayTDiff>0 && fImplantMap[temp]->IsGood()){ // decay must came later than imp;
        double ionDecayTDiff = fImplantMap[temp]->fPIN1T - T; 
        if(ionDecayTDiff <dt) {
          dt = ionDecayTDiff;
          best = temp;
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
    //FillHistogram("pid_cur", 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    //if(current->Stopped()){
    //  FillHistogram("pidnosssd_cur", 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    //}


    //TIter iter(gCuts);
    //while(TCutG *cut = (TCutG*)iter.Next()) {
    //  if(!cut->IsInside(current->fI2S,current->fPIN1E)) continue;    
    //  FillHistogram("pid",cut, 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    //  if(current->Stopped()) 
    //    FillHistogram("pidnosssd",cut, 2e3,0,2e4,current->fI2S, 8e3,0,8e3,current->fPIN1E);
    //  for(auto &it1 : *vdec) {

    //     FillHistogram("summary",cut, 300,0,300,it1.fDSSDFront[0].GetNumber(), 8000,0,4000,it1.fDSSDFront[0].GetEnergy());
    //     FillHistogram("summary_raw",cut, 300,0,300,it1.fDSSDFront[0].GetNumber(), 8000,0,64000,it1.fDSSDFront[0].GetCharge());
    //     FillHistogram("summary",cut, 300,0,300,it1.fDSSDBack[0].GetNumber(), 8000,0,4000,it1.fDSSDBack[0].GetEnergy());
    //     FillHistogram("summary_raw",cut, 300,0,300,it1.fDSSDBack[0].GetNumber(), 8000,0,64000,it1.fDSSDBack[0].GetCharge());


    //    for(auto &it2 : it1.fGe) { 
    //      double dt = it1.GetTimestamp() - it2.GetTimestamp();

    //      FillHistogram("timedif_decay_gamma",cut,200,-1000,1000,dt);
    //      FillHistogram("timedif2_decay_gamma",cut,200,-1000,1000,dt,4000,0,4000,it2.GetEnergy());

    //      FillHistogram("summary",cut, 300,0,300,it2.GetNumber(), 8000,0,4000,it2.GetEnergy());
    //      FillHistogram("summary_raw",cut, 300,0,300,it2.GetNumber(), 8000,0,64000,it2.GetCharge());
 
    //      //double deltatime = fabs(it1.GetTimestamp() - it2.GetTimestamp());

    //      //if(deltatime>40 && deltatime<300){ // Fill gamma spectrum
    //      if(dt>-200 && dt<300){ // Fill gamma spectrum
    //        FillHistogram("summary_prompt",cut, 300,0,300,it2.GetNumber(), 8000,0,4000,it2.GetEnergy());
    //        FillHistogram("summary_prompt_raw",cut, 300,0,300,it2.GetNumber(), 8000,0,64000,it2.GetCharge());
    //      }
    //    }
    //  }
    //}

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
