#include "DetHit.h"
#include<iostream>
#include <ddaschannel.h>

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
double BCSEvent::I2S() const  { for(auto &it : fHits) { if(it.GetNumber()==177) return it.GetCharge(); } return -1; }


int BCSEvent::HGFSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),0,39))    count++; } return count; }
int BCSEvent::HGBSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),80,119))  count++; } return count; }
int BCSEvent::LGFSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),40,79))   count++; } return count; }
int BCSEvent::LGBSize() const {  int count=0; for(auto &it : fHits) { if(Range(it.GetNumber(),120,159)) count++; } return count; }




