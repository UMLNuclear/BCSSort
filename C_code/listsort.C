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
  std::vector<DetHit> FL;
  std::vector<DetHit> BL;
  DetHit pin1_i2s;
  DetHit pin1;
  DetHit i2s_i2n;
  double pin1_charge = -1;
  double pin1_cout = 0;
  for(size_t y=0;y<vec.size();y++){
    if(vec[y].GetNumber()==181){
      pin1_cout++;
      FillHistogram("PIN1E_in",10e3,0,10e3,vec[y].GetCharge());
      if(vec[y].GetCharge()>100) FillHistogram("PIN1E2_in",10e3,0,10e3,vec[y].GetCharge());
    }
    FillHistogram("PIN1_Count",10,0,10,pin1_cout);
    switch(vec[y].GetNumber()){
      case 40 ... 79:  //FL
        FL.push_back(vec[y]);
        break;

      case 120 ... 159:  //BL
        BL.push_back(vec[y]);
        break;

      case 160 ... 175:  //SSSD
        if(vec[y].GetCharge()>100){
          sssd_flag++;
        }
        break;

      case 177: //pin1_i2s
        pin1_i2s = vec[y];
        break;

      case 180: //i2s_i2n
        i2s_i2n = vec[y];
        break;

      case 181: //pin1
        pin1 = vec[y];
        pin1_charge = vec[y].GetCharge();
        break;

      default:
        break;
    }
  }
  if(pin1_charge>0){
    FillHistogram("PIN1E",10e3,0,10e3,pin1_charge);
    if(pin1_charge>100) FillHistogram("PIN1E2",10e3,0,10e3,pin1_charge);
    if(sssd_flag==0){
      if(FL.size()>0 && BL.size()>0){
        //FillHistogram("pid_isgood_nosssd",2e3,0,2e4,pin1_i2s.GetCharge(), 10e3,0,10e3,pin1.GetCharge());
        FillHistogram("PIN1E_isgood_nosssd", 10e3,0,10e3,pin1.GetCharge());
        if(pin1.GetCharge()>100) FillHistogram("PIN1E_isgood_nosssd2", 10e3,0,10e3,pin1.GetCharge());
      }
      FillHistogram("PIN1E_nosssd",10e3,0,10e3,pin1.GetCharge());
      if(pin1.GetCharge()>100) FillHistogram("PIN1E_nosss2",10e3,0,10e3,pin1.GetCharge());
    }  
  }


}




void listsort(){

  int faddress;
  int fnumber;
  double ftimestamp;
  double fcharge;
  gChain->SetBranchAddress("address",   &faddress);
  gChain->SetBranchAddress("number",    &fnumber);
  gChain->SetBranchAddress("timestamp", &ftimestamp);
  gChain->SetBranchAddress("charge",    &fcharge);
  TChannel::ReadDetMapFile();

  double starttime = 0;
  double buildtime = BUILDTIME;
  std::vector<DetHit> vec_hit;

  long n = gChain->GetEntries();
  //n = 1e6;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    DetHit *fhit = new DetHit;
    fhit->SetAddress(faddress);
    fhit->SetNumber(fnumber);
    fhit->SetTimestamp(ftimestamp);
    fhit->SetCharge(fcharge);
    if(fhit->GetNumber()==181){
      FillHistogram("PIN1E_sin",10e3,0,10e3,fhit->GetCharge());
      if(fhit->GetCharge()>100) FillHistogram("PIN1E2_sin",10e3,0,10e3,fhit->GetCharge());
    }
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
  SaveHistograms("listoutput1031.root");
  return;
}
