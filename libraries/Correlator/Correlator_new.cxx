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

void Correlator::AddImplant(Implant *fimp){
  
  pixel pix = fimp->GetPixel();
  if(pix.first<0 || pix.first>39) return;
  if(pix.second<0 || pix.second>39) return;
  if(fBlockMap[pix]>0){ // pixel was implanted
    fBlockMap[pix] = fimp->fPIN1T;
    FlushPixel(pix,false);
  }
  fimp->Copy(*fImplantMap[pix]);
  return;
}

void AddDecay(Decay *fdec){
  pixel pix = fimp->GetPixel();
  if(pix.first<0 || pix.first>39) return;
  if(pix.second<0 || pix.second>39) return;
  Match(fdec);
}


void Match(Decay *fdec, bool flagT = false){
  double dt = DBL_MAX;
  pixel best = std::make_pair(-1,-1);
  pixel pix = fdec->GetPixel();
  for(int i=pix.first-1;i<=pix.first+1;i++){
    if(i<0 || i>39) continue;
    for(int j=pix.first-1;j<=pix.first+1;j++){
      if(j<0 || j>39) continue;
      pixel temp = std::make_pair(i,j);
      if(fBlockMap[temp]>0 && fBlockMap[temp]<fdec->DSSDhiT()){ //pixel was implanted
        if(!flagT){
          fDecayMap[temp]->push_back(*fDecay);
        }else{
          if(fabs(fImplantMap[temp]->fPIN1T - fdec->DSSDhiT())<dt){
            dt = fImplantMap[temp]->fPIN1T - fdec->DSSDhiT();
            best = temp;
          }
        }
      }
    }
  }
  if(flagT){
    if(best.first<0 || best.second<0) continue;
    if(best.first>39 || best.second>39) continue;
    fDecayMap[best]->push_back(*fDecay);
  }  
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
