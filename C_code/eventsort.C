#include<map>


#include<TFile.h>
#include<TCutG.h>
#include<TH2D.h>
#include<TChain.h>
#include<TKey.h>
#include<TList.h>
#include<TCutG.h>

#include<Implant.h>
#include<util.h>
#include<TChannel.h>
#include<DetHit.h>
#include<TOFCorrection.h>


std::map<pixel, double> fChargeMap; //pix = <FL#, BL#>, DetHit = pin1
std::map<pixel, double> fTimeMap; //pix = <FL#, BL#>, DetHit = pin1



std::pair<int, int> GetPixel(std::vector<DetHit> LGF, std::vector<DetHit> LGB){
  int indexf = (-LGF[0].GetNumber()+80)-1;
  int indexb = (-LGB[0].GetNumber()+160)-1;
  double tempE = LGF[0].GetEnergy();
  for(int z=1;z<LGF.size();z++){
    if(tempE<LGF[z].GetEnergy()){
      tempE = LGF[z].GetEnergy();
      indexf = (-LGF[z].GetNumber()+80)-1;   // from 0~39: FL1->0; FL40->39
    }
  }
  tempE = LGB[0].GetEnergy();
  for(int z=1;z<LGB.size();z++){
    if(tempE<LGB[z].GetEnergy()){
      tempE = LGB[z].GetEnergy();
      indexb = (-LGB[z].GetNumber()+160)-1;
    }
  }
  return std::make_pair(indexf,indexb);
}


void eventsort(){

  for(int m=0;m<40;m++){
    for(int n=n;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fChargeMap[pix] = -1;
      fTimeMap[pix] = -1;
    }
  }



  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    DetHit pin1;
    for(auto &it:fevent->fHits){
      if(fevent->Pin1E()>0 && fevent->LGFSize()==0 && fevent->LGBSize()==0){
        FillHistogram("summary",16e3,0,8e4,it.GetCharge(), 300,0,300,it.GetNumber());
        FillHistogram("summary_cal",16e3,0,8e4,it.GetEnergy(), 300,0,300,it.GetNumber());
      }
      if(it.GetNumber()==181){
        pin1 = it;
        if(pin1.GetCharge()>0){

          if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
            FillHistogram("PIN1E_isgood_sssd", 10e3,0,10e3,pin1.GetCharge());
            if(pin1.GetCharge()>100){
              FillHistogram("PIN1E_isgood_sssd2", 10e3,0,10e3,pin1.GetCharge());
              FillHistogram("PIN1Esize_isgood_sssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
              FillHistogram("ESize_isgood_sssd", 100,0,100,fevent->Size());
              for(auto &it1:fevent->fHits){
                FillHistogram("Map_isgood_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
              }
            }

            std::pair<int, int> pix = GetPixel(fevent->LGF(), fevent->LGB());
            if(fTimeMap[pix]<0){
              fChargeMap[pix] = pin1.GetCharge();
              fTimeMap[pix] = pin1.GetTimestamp();
            }else{
              double dtime = fabs(fTimeMap[pix] - pin1.GetTimestamp());
              double pin1c = fChargeMap[pix];
              if(dtime>500e6){
                FillHistogram("PIN1E_block_sssd5",10e3,0,10e3,pin1c);
                if(pin1c>100){
                  FillHistogram("PIN1E2_block_sssd5",10e3,0,10e3,pin1c);
                  FillHistogram("PIN1Esize_block_sssd5", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
                  FillHistogram("ESize_block_sssd5", 100,0,100,fevent->Size());
                  for(auto &it1:fevent->fHits){
                    FillHistogram("Map_block_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                  }
                }                
              }  
              fChargeMap[pix] = pin1.GetCharge();
              fTimeMap[pix] = pin1.GetTimestamp();
            }
          }
          if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
            FillHistogram("PIN1E_ForB_sssd", 10e3,0,10e3,pin1.GetCharge());
            if(pin1.GetCharge()>100){
              FillHistogram("PIN1E_ForB_sssd2", 10e3,0,10e3,pin1.GetCharge());
              FillHistogram("PIN1Esize_ForB_sssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
              FillHistogram("ESize_ForB_sssd", 100,0,100,fevent->Size());
              for(auto &it1:fevent->fHits){
                if(pin1.GetCharge()>5000 && pin1.GetCharge()<6000){
                  FillHistogram("Map_ForB_Ne_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                  FillHistogram("Sum_ForB_Ne_sssd",300,0,300,it1.GetNumber(),30e2,0,30e3,it1.GetEnergy());
                }
                if(pin1.GetCharge()>6000 && pin1.GetCharge()<7300){
                  FillHistogram("Map_ForB_Na_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                  FillHistogram("Sum_ForB_Na_sssd",300,0,300,it1.GetNumber(),30e2,0,30e3,it1.GetEnergy());
                }
              }
            }
          } 
          if(fevent->LGFSize()==0 && fevent->LGBSize()==0){
            FillHistogram("PIN1E_noLG_sssd", 10e3,0,10e3,pin1.GetCharge());
            if(pin1.GetCharge()>100){

              FillHistogram("PIN1E_noLG_sssd2", 10e3,0,10e3,pin1.GetCharge());
              FillHistogram("PIN1Esize_noLG_sssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
              FillHistogram("ESize_noLG_sssd", 100,0,100,fevent->Size());
              if(pin1.GetCharge()>5000 && pin1.GetCharge()<6000){
                if(fevent->HGFSize()>0) FillHistogram("HGFSize_Ne",50,0,50,fevent->HGFSize());
                if(fevent->HGFSize()==0 && fevent->HGBSize()==0 && fevent->SSSDSize()>0) 
                  FillHistogram("NeSSSD_NoDSSD", 20,0,20,fevent->SSSDSize());
              }
              if(pin1.GetCharge()>6000 && pin1.GetCharge()<7300){
                if(fevent->HGFSize()>0)FillHistogram("HGFSize_Na",50,0,50,fevent->HGFSize());
                if(fevent->HGFSize()==0 && fevent->HGBSize()==0 && fevent->SSSDSize()>0) 
                  FillHistogram("NaSSSD_NoDSSD", 20,0,20,fevent->SSSDSize());
              }
              for(auto &it1:fevent->fHits){
                if(pin1.GetCharge()>5000 && pin1.GetCharge()<6000){
                  FillHistogram("Map_noLG_Ne_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                  FillHistogram("Sum_noLG_Ne_sssd",300,0,300,it1.GetNumber(),30e2,0,30e3,it1.GetEnergy());
                }
                if(pin1.GetCharge()>6000 && pin1.GetCharge()<7300){
                  FillHistogram("Map_noLG_Na_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                  FillHistogram("Sum_noLG_Na_sssd",300,0,300,it1.GetNumber(),30e2,0,30e3,it1.GetEnergy());
                }
              }
            }
          } 
          FillHistogram("PIN1E_sssd",10e3,0,10e3,pin1.GetCharge());
          if(pin1.GetCharge()>100){
            FillHistogram("PIN1E2_sssd",10e3,0,10e3,pin1.GetCharge());
            FillHistogram("PIN1Esize_sssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
            FillHistogram("ESize_sssd", 100,0,100,fevent->Size());
            for(auto &it1:fevent->fHits){
              FillHistogram("Map_sssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
            }
          }

          //===================================== NO SSSD ====================================================//
          if(fevent->SSSDSize()==0){
            if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
              FillHistogram("PIN1E_isgood_nosssd", 10e3,0,10e3,pin1.GetCharge());
              if(pin1.GetCharge()>100){
                FillHistogram("PIN1E_isgood_nosssd2", 10e3,0,10e3,pin1.GetCharge());
                FillHistogram("PIN1Esize_isgood_nosssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
                FillHistogram("ESize_isgood_nosssd", 100,0,100,fevent->Size());
                for(auto &it1:fevent->fHits){
                  FillHistogram("Map_isgood_nosssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                }
              }


              std::pair<int, int> pix = GetPixel(fevent->LGF(), fevent->LGB());
              if(fTimeMap[pix]<0){
                fChargeMap[pix] = pin1.GetCharge();
                fTimeMap[pix] = pin1.GetTimestamp();
              }else{
                double dtime = fabs(fTimeMap[pix] - pin1.GetTimestamp());
                double pin1c = fChargeMap[pix];
                if(dtime>500e6){
                  FillHistogram("PIN1E_block_nosssd5",10e3,0,10e3,pin1c);
                  if(pin1c>100){
                    FillHistogram("PIN1E2_block_nosssd5",10e3,0,10e3,pin1c);
                    FillHistogram("PIN1Esize_block_nosssd5", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
                    FillHistogram("ESize_block_nosssd5", 100,0,100,fevent->Size());
                    for(auto &it1:fevent->fHits){
                      FillHistogram("Map_block_nosssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                    }
                  }                
                }  
                fChargeMap[pix] = pin1.GetCharge();
                fTimeMap[pix] = pin1.GetTimestamp();
              }



            }
            if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
              FillHistogram("PIN1E_ForB_nosssd", 10e3,0,10e3,pin1.GetCharge());
              if(pin1.GetCharge()>100){
                FillHistogram("PIN1E_ForB_nosssd2", 10e3,0,10e3,pin1.GetCharge());
                FillHistogram("PIN1Esize_ForB_nosssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
                FillHistogram("ESize_ForB_nosssd", 100,0,100,fevent->Size());
                for(auto &it1:fevent->fHits){
                  FillHistogram("Map_ForB_nosssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
                }
              }
            } 
            FillHistogram("PIN1E_nosssd",10e3,0,10e3,pin1.GetCharge());
            if(pin1.GetCharge()>100){
              FillHistogram("PIN1E2_nosssd",10e3,0,10e3,pin1.GetCharge());
              FillHistogram("PIN1Esize_nosssd", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
              FillHistogram("ESize_nosssd", 100,0,100,fevent->Size());
              for(auto &it1:fevent->fHits){
                FillHistogram("Map_nosssd",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
              }
            }
          } 

          //===============================================================================================// 
          FillHistogram("PIN1E",10e3,0,10e3,pin1.GetCharge());
          if(pin1.GetCharge()>100){
            FillHistogram("PIN1E2",10e3,0,10e3,pin1.GetCharge());
            FillHistogram("PIN1Esize", 10e3,0,10e3,pin1.GetCharge(), 100,0,100,fevent->Size());
            FillHistogram("ESize", 100,0,100,fevent->Size());
            for(auto &it1:fevent->fHits){
              FillHistogram("Map",300,0,300,it1.GetNumber(),100,0,100,fevent->Size());
            }
          }
        }
      }
    }   


    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms("eventoutput1031.root");
  return;
}
