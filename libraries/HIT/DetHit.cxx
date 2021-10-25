#include "DetHit.h"
#include<iostream>
#include <ddaschannel.h>
#include <TDirectory.h>
#include <globals.h>
#include <numeric>

DetHit::DetHit() {  }//create, constructor

DetHit::DetHit(const DetHit &rhs) {
    timestamp = rhs.GetTimestamp();
    charge    = rhs.GetCharge();
    number    = rhs.GetNumber();
    address   = rhs.GetAddress();
}


DetHit::DetHit(ddaschannel *chan) {
    timestamp = chan->GetCoarseTime();
    //int crate  = chan->GetCrateID();
    //int slot   = chan->GetSlotID(); 
    //int nch    = chan->GetChannelID(); 
    charge  = chan->Energy();
    number  = chan->GetNumber();
    address = chan->GetAddress();
}

DetHit::~DetHit() {  }//delete, destructor




void DetHit::print(){
    std::cout<<number <<"\t"<< charge <<"\t"<< timestamp<<std::endl;
}



double DetHit::GetEnergy() const {
    TChannel *c = TChannel::Get(address);
    //cout << c->GetNumber() <<"\t" 
    if(!c) return -1;
    //if(abs((GetNumber()-240)<32)) {
    //  std::cout << GetNumber() <<"\t" << GetCharge() << "\t" << GetEnergy() << std::endl;
    //}

    return c->CalEnergy(charge); 
}


BCSEvent::BCSEvent() { }

BCSEvent::~BCSEvent() { }

double BCSEvent::Pin1E() const  { for(auto &it : fHits) { if(it.GetNumber()==181) return it.GetCharge();    } return -1; }
double BCSEvent::Pin1T() const  { for(auto &it : fHits) { if(it.GetNumber()==181) return it.GetTimestamp(); } return -1; }
double BCSEvent::Pin2E() const  { for(auto &it : fHits) { if(it.GetNumber()==182) return it.GetCharge();    } return -1; }
double BCSEvent::Pin2T() const  { for(auto &it : fHits) { if(it.GetNumber()==182) return it.GetTimestamp(); } return -1; }
double BCSEvent::Pin3E() const  { for(auto &it : fHits) { if(it.GetNumber()==183) return it.GetCharge();    } return -1; }
double BCSEvent::Pin3T() const  { for(auto &it : fHits) { if(it.GetNumber()==183) return it.GetTimestamp(); } return -1; }
double BCSEvent::I2S() const  { for(auto &it : fHits) { if(it.GetNumber()==177) return it.GetCharge(); } return -1; }
double BCSEvent::PIN1_I2N() const  { for(auto &it : fHits) { if(it.GetNumber()==176) return it.GetCharge(); } return -1; }
double BCSEvent::I2N_I2S() const  { for(auto &it : fHits) { if(it.GetNumber()==180) return it.GetCharge(); } return -1; }

int BCSEvent::HGFSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),0,39))    count++; } return count; }
int BCSEvent::HGBSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),80,119))  count++; } return count; }
int BCSEvent::LGFSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),40,79))   count++; } return count; }
int BCSEvent::LGBSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),120,159)) count++; } return count; }
int BCSEvent::SSSDSize() const {  
  int count=0; 
  for(auto &it : fHits) { 
    if(Range(it.GetNumber(),160,175)){
      if(it.GetCharge()>100) count++; 
    }  
  } 
  return count; 
}




std::vector<DetHit> BCSEvent::LGF() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),40,79)) vec.push_back(it);
  } 
  return vec;
}

std::vector<DetHit> BCSEvent::LGB() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),120,159)) vec.push_back(it);
  } 
  return vec;
}

int BCSEvent::LGPixel() const {

  std::vector<DetHit> lgf = LGF();
  std::vector<DetHit> lgb = LGB();
  if(lgf.size()!=1 || lgb.size()!=1) return -1;
  return lgf.at(0).GetNumber()*1000 + lgb.at(0).GetNumber();

}

std::vector<DetHit> BCSEvent::HGF() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),0,39)) vec.push_back(it);
  } 
  return vec;
}

std::vector<DetHit> BCSEvent::HGB() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),80,119)) vec.push_back(it);
  } 
  return vec;
}

//============================ return HG Pixel =========================//
std::pair<int,int> BCSEvent::HGPixel() const {

  int x = -1;
  int y = -1;
  if(HGFSize()<2 || HGBSize()<2){
    return std::make_pair(x,y);
  }
  if(HGFSize()>20 || HGBSize()>20){
    return std::make_pair(x,y);
  }

  std::vector<DetHit> hgf = HGF();
  std::vector<std::pair<double, int>> fvec; // <charge, number>;
  for(size_t i=0;i<HGFSize();i++){
    double charge = hgf[i].GetCharge();
    int num = hgf[i].GetNumber();
    fvec.push_back(std::make_pair(charge, num));
  }
  std::sort(fvec.begin(),fvec.end());
  std::reverse(fvec.begin(),fvec.end());
  std::vector<int> fnumvec;
  bool fflag = false;
  for(size_t i=0;i<fvec.size();i++){
    //if(fvec[i].first<(0.9*fvec[0].first)) break; // charge < 90% * charge_max => bad;
    for(size_t j=i+1;j<fvec.size();j++){
      if(fabs(fvec[i].second - fvec[j].second)<2){//max energy hit has neighbors
        fflag = true;
        double cdif = fvec[i].first - fvec[j].first;
        if(cdif<(0.1*fvec[i].first)){ //its neighbor's charge is close to it
          fnumvec.push_back(fvec[j].second);
        }
      }
    }
    if(fflag){
      fnumvec.insert(fnumvec.begin(),fvec[i].second);
      break;
    } 
  }
  if(fnumvec.size()==0) x=-1;// good channel must have neighbor;
  else x = std::accumulate(fnumvec.begin(),fnumvec.end(),0)/fnumvec.size();
  
  std::vector<DetHit> hgb = HGB();
  std::vector<std::pair<double, int>> bvec; // <charge, number>;
  for(size_t i=0;i<HGBSize();i++){
    double charge = hgb[i].GetCharge();
    int num = hgb[i].GetNumber();
    bvec.push_back(std::make_pair(charge, num));
  }
  std::sort(bvec.begin(),bvec.end());
  std::reverse(bvec.begin(),bvec.end());
  std::vector<int> bnumvec;
  bool bflag = false;
  for(size_t i=0;i<bvec.size();i++){
    //if(bvec[i].first<(0.9*bvec[0].first)) break; // charge < 10% * charge_max => bad;
    for(size_t j=i+1;j<bvec.size();j++){
      if(fabs(bvec[i].second - bvec[j].second)<2){//max energy hit has neighbors
        bflag = true;
        double cdif = bvec[i].first - bvec[j].first;
        if(cdif<(0.1*bvec[i].first)){ //its neighbor's charge is close to it
          bnumvec.push_back(bvec[j].second);
        }
      }
    }
    if(bflag){
      bnumvec.insert(bnumvec.begin(),bvec[i].second);
      break;
    } 
  }
  if(bnumvec.size()==0) y=-1;
  else y = (std::accumulate(bnumvec.begin(),bnumvec.end(),0)/bnumvec.size())-80;
  
  //for(size_t j=0;j<bvec.size();j++){
  //  printf("x[%02i]: %f\n",bvec[j].second-80, bvec[j].first);
  //}
  //for(size_t i=0;i<fvec.size();i++){
  //  printf("y[%02i] : %f\n",fvec[i].second, fvec[i].first);
  //}

  return std::make_pair(y,x);
}
//===================================================================//

TH2D *BCSEvent::DrawHG(Option_t *opt) const {

  TH2D *hist = 0; 
  gDirectory->FindObjectAny("hitpad");
  if(!hist) {
    hist = new TH2D("hitpad","hitpad",80,0,80, 40,0,40);
  } else {
    hist->Reset();
  }
  bool print=false;
  TString sopt(opt);
  sopt.ToLower();
  if(sopt.Contains("print")) print = true;
  for(size_t m=0;m<HGFSize();m++){
    for(size_t n=0;n<HGBSize();n++){
      int x,y, fc,bc;
       x = HGB().at(n).GetNumber()-80;
       y = HGF().at(m).GetNumber();
       fc = HGF().at(m).GetCharge();
       bc = HGB().at(n).GetCharge();
       if(print) printf("\t x[%02i]\ty[%02i] : %i\t%i\n",x,y,bc,fc); 
       hist->Fill(x,y,fc);
       hist->Fill(x+40,y,bc);
    }
  }
  
  return hist;
}


void BCSEvent::PrintNumber(){
  int clgf, clgb, chgf, chgb = 0; 
  int csssd, clabr, chpge, cig = 0;
  std::vector<std::pair<int, double>> vec;
  printf("Event Size = %d\n\n", Size());
  for(auto &it:fHits){  
    int num = it.GetNumber();
    if(num>=0 && num<=39) chgf++;
    if(num>=40 && num<=79) clgf++;
    if(num>=80 && num<=119) chgb++;
    if(num>=120 && num<=159) clgb++;
    if(num>=160 && num<=175) csssd++;
    if(num>=192 && num<=207) clabr++;
    if(num>=208 && num<=271) chpge++;
    if(num>=272 && num<=287) cig++;
    vec.push_back(std::make_pair(it.GetNumber(), it.GetEnergy()));
  }
  std::sort(vec.begin(), vec.end());
  for(size_t i=0;i<vec.size();i++){
    printf("Detector # = %i\t\t", vec[i].first);
    printf("Energy  = %f\n", vec[i].second);
  }
  printf("\n");
  printf("LGF.size = %i\n",clgf);
  printf("LGB.size = %i\n",clgb);
  printf("HGF.size = %i\n",chgf);
  printf("HGB.size = %i\n",chgb);
  printf("SSSD.size = %i\n",csssd);
  printf("HPGe.size = %i\n",chpge);
  printf("LaBr.size = %i\n",clabr);
  printf("SSSD_Ig.size = %i\n",cig);
  printf("====================================\n");
}



