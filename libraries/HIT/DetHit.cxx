#include "DetHit.h"
#include<iostream>
#include <ddaschannel.h>
#include <TDirectory.h>
#include <globals.h>
#include <numeric>

DetHit::DetHit() { energy = sqrt(-1);  }//create, constructor

DetHit::DetHit(const DetHit &rhs) {
    timestamp = rhs.GetTimestamp();
    charge    = rhs.GetCharge();
    number    = rhs.GetNumber();
    address   = rhs.GetAddress();
    energy    = rhs.GetEnergy();
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



double DetHit::GetEnergy(bool recal) const {
    TChannel *c = TChannel::Get(address);
    //cout << c->GetNumber() <<"\t" 
    if(!c) return -1;
    //if(abs((GetNumber()-240)<32)) {
    //  std::cout << GetNumber() <<"\t" << GetCharge() << "\t" << GetEnergy() << std::endl;
    //}
    if(energy!=energy || recal) 
      energy = c->CalEnergy(charge); 

    return energy; 
}



void DetHit::Clear(){
    timestamp = -1;
    charge    = -1;
    number    = -1;
    address   = -1;
    energy    = sqrt(-1);
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
int BCSEvent::HPGeSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),208,271)) count++; } return count; }
int BCSEvent::SSSDSize() const {  
  int count=0; 
  for(auto &it : fHits) { 
    if(Range(it.GetNumber(),160,175)){
      //if(it.GetCharge()>100) count++;
      count++; 
    }  
  } 
  return count; 
}
int BCSEvent::SSSDLGSize() const {  
  int count=0; 
  for(auto &it : fHits) { 
    if(Range(it.GetNumber(),272,287)){
      //if(it.GetCharge()>100) 
      count++; 
    }  
  } 
  return count; 
}


DetHit BCSEvent::LGFMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:LGF()){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

DetHit BCSEvent::LGBMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:LGB()){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

DetHit BCSEvent::HGFMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:HGF()){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

DetHit BCSEvent::HGBMax() const {
  DetHit fhit;
  double Emax = 0;
  for(auto &it:HGB()){
    if(it.GetEnergy()>Emax){
      fhit = it;
      Emax = it.GetEnergy();
    }
  }
  return fhit;
}

double BCSEvent::DSSDloT() const {
  double t = -1;
  if(LGFSize()>0){
    t = LGFMax().GetTimestamp();    
  }else if(LGBSize()>0){
    t = LGBMax().GetTimestamp();
  }
  return t;
}

double BCSEvent::DSSDhiT() const {
  double t = -1;
  if(HGFSize()>0){
    t = HGFMax().GetTimestamp();    
  }else if(HGBSize()>0){
    t = HGBMax().GetTimestamp();
  }
  return t;
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

std::vector<DetHit> BCSEvent::HPGe() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),208,271)) vec.push_back(it);
  } 
  return vec;
}

std::vector<DetHit> BCSEvent::LaBr() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),192,207)) vec.push_back(it);
  } 
  return vec;
}


std::vector<DetHit> BCSEvent::SSSD() const {
  std::vector<DetHit> vec;
  for(auto &it :fHits) {
    if(Range(it.GetNumber(),160,175)) vec.push_back(it);
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

TH2D *BCSEvent::DrawLG(Option_t *opt) const {

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
  for(size_t m=0;m<LGFSize();m++){
    for(size_t n=0;n<LGBSize();n++){
      int x,y, fc,bc;
       x = -LGB().at(n).GetNumber()+159;
       y = -LGF().at(m).GetNumber()+79;
       fc = LGF().at(m).GetEnergy();
       bc = LGB().at(n).GetEnergy();
       if(print) printf("\t x[%02i]\ty[%02i] : %i\t%i\n",x,y,bc,fc); 
       hist->Fill(x,y,fc);
       hist->Fill(x+40,y,bc);
    }
  }
  
  return hist;
}

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


TH2D *BCSEvent::DrawSSSD(Option_t *opt) const {
  
  TH2D *hist = 0;
  gDirectory->FindObjectAny("sum");
  if(!hist){
    hist = new TH2D("sum","sum",16,0,16, 8e3,0,8e3);
  }else{
    hist->Reset();
  }
  bool print = false;
  TString sopt(opt);
  sopt.ToLower();
  if(sopt.Contains("print")) print = true;
  for(size_t m=0;m<SSSDSize();m++){
    int x = SSSD().at(m).GetNumber()-160;
    double charge = SSSD().at(m).GetCharge();
    if(print) printf("\t channel = %02i\tcharge = %f\n",x,charge);
    hist->Fill(x,charge);
  }
}


//===================================================//
void BCSEvent::Print(){
  std::cout<<"number"<<"\t"<<"energy"<<"\t"<<"time"<<std::endl;
  for(auto &it:fHits){
    std::cout<<it.GetNumber()<<"\t"<<it.GetEnergy()<<"\t"<<it.GetTimestamp()<<std::endl;
  }
}











