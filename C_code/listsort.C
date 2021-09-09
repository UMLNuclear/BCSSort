#include<map>


#include<TFile.h>
#include<TCutG.h>
#include<TH2D.h>
#include<TChain.h>
#include<TKey.h>
#include<TList.h>

#include<Implant.h>
#include<util.h>
#include<TChannel.h>
#include<DetHit.h>



void Process(std::vector<DetHit> vec){
  int sssd_flag = 0;
  DetHit pin1_i2s;
  DetHit pin1;
  DetHit i2s_i2n;
  for(size_t x=0;x<vec.size();x++){
    switch(vec[x].GetNumber()){
      case 160 ... 175:  //SSSD
        if(vec[x].GetCharge()>0){
          sssd_flag++;
        }
        break;
      
      case 177: //pin1_i2s
        pin1_i2s = vec[x];
        break;
      
      case 180: //i2s_i2n
        i2s_i2n = vec[x];
        break;
      
      case 181: //pin1
        pin1 = vec[x];
        break;

        default:
        break;
    }
    if(pin1_i2s.GetEnergy()>0 && pin1.GetCharge()>0){
      if(!sssd_flag){
        FillHistogram("pid_nosssd",3e3,0,3e4,pin1_i2s.GetCharge(), 8e3,0,2e4,pin1.GetCharge());
      }  
      FillHistogram("pid",3e3,0,3e4,pin1_i2s.GetCharge(), 8e3,0,2e4,pin1.GetCharge());
    }
    if(pin1_i2s.GetEnergy()>0 && i2s_i2n.GetCharge()>0){
      if(!sssd_flag){
        FillHistogram("momentum_nosssd",3e3,0,3e4,pin1_i2s.GetCharge(), 3e3,0,3e4,i2s_i2n.GetCharge());
      }
      FillHistogram("momentum",3e3,0,3e4,pin1_i2s.GetCharge(), 3e3,0,3e4,i2s_i2n.GetCharge());
    }
    
  }



}




void listsort(){

  //DetHit *fhit = new DetHit;
  DetHit *fhit;
  gChain->SetBranchAddress("DetHit", &fhit);
  TChannel::ReadDetMapFile();

  double starttime = 0;
  double buildtime = BUILDTIME;
  std::vector<DetHit> vec_hit;
  
  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if((fhit->GetTimestamp()-starttime) > buildtime){
      starttime = fhit->GetTimestamp();
      Process(vec_hit);
      vec_hit.clear();
    }
    vec_hit.push_back(*fhit);
    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms("listoutput1100.root");
  return;
}
