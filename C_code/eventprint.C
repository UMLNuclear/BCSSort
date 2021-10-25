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
#include<globals.h>



void eventprint(){

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;
  long c1=0;
  long c2=0;
  long c3=0;
  long c4=0;
  long c5=0;
  long c6=0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    //if(fevent->Pin1E()>5000 && fevent->Pin1E()<6000){ 
    if(fevent->Pin1E()>100){ 
      if(fevent->LGFSize()==0 && fevent->LGBSize()==0){
        c1++;
        pixel pix = fevent->HGPixel();
        if(pix.first>0 && pix.second>0){
          c2++;
        }
      }
    }
    if(fevent->Pin1E()>5000 && fevent->Pin1E()<6000){ 
      if(fevent->LGFSize()==0 && fevent->LGBSize()==0){
        c3++;
        pixel pix = fevent->HGPixel();
        if(pix.first>0 && pix.second>0){
          c4++;
        }
      }
    }
    if(fevent->Pin1E()>6000 && fevent->Pin1E()<7300){ 
      if(fevent->LGFSize()==0 && fevent->LGBSize()==0){
        c5++;
        pixel pix = fevent->HGPixel();
        if(pix.first>0 && pix.second>0){
          c6++;
        }
      }
    }
  }
  printf("NoLG counts = %lu\n", c1); 
  printf("NoLG+GoodHG counts = %lu\n", c2); 
  printf("Ne NoLG counts = %lu\n", c3); 
  printf("Ne NoLG+GoodHG counts = %lu\n", c4); 
  printf("Na NoLG counts = %lu\n", c5); 
  printf("Na NoLG+GoodHG counts = %lu\n", c6); 
  //SaveHistograms("junk.root");
  return;

}

