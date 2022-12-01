

//#include "global.h"
//#include "DetHit.h"
//#include "TOFCOrrection.h"

typedef std::pair<int, int> pixel;

std::vector<double> tofpar; 
std::vector<double> tofpar1;

TCutG *ban;
TCutG *cutstrip;
TCutG *blob[6];
TH2D *h1 = 0;
TH2D *h2[6];
TH1D *h3[6];
TH1D *h4[6];

std::map<pixel, BCSEvent*>fImpMap;
std::map<pixel, double> fBLMap;
std::map<pixel, int> fDecCounter;

int count0 = 0; // counts out of veto;
int count1 = 0; // good implant will be written in TTree;
int count2 = 0; // black list: back to back;
int count3 = 0; // black list: pixel is still cooling down;  

bool veto(DetHit lgf, std::vector<DetHit> sssd){

  int    lgfn = lgf.GetStrip()-1;
  double lgfe = lgf.GetEnergy();

  for(int i=0;i<sssd.size();i++){
    int sssdn  = sssd[i].GetStrip()-1;
    long sssdc = sssd[i].GetCharge();

    if(cutstrip->IsInside(sssdn, lgfn)){
      if(ban->IsInside(sssdc, lgfe)){
        return true;
      }
    }
  }
  return false;
}


void InitMap(){ 
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new BCSEvent;
      fBLMap[pix] = -1;
      fDecCounter[pix] = 0;
    }
  }
}

void CleanMap(double pin1t){
  for(auto &it:fImpMap){
    if(it.second->Pin1T()<0) continue; // this pixel is not implanted;
    double dt1 = pin1t - it.second->Pin1T();
    if(dt1>1e9){
      //TODO::it.second->Fill TTree;
      double pin1e = it.second->Pin1E();
      double tof   = it.second->I2S();
      tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
      h1->Fill(tof, pin1e);
      count1++;
      for(int m=0;m<6;m++){
        if(blob[m]->IsInside(tof,pin1e)){
          h4[m]->Fill(fDecCounter[it.first]);
        }
      }
      it.second->Clear();
      fDecCounter[it.first]=0;
    }
  }

  for(auto &it:fBLMap){
    if(it.second<0) continue; // this pixel is not blocked;
    double dt2 = pin1t - it.second;
    if(dt2>1e9){ 
      it.second = -1; //unblock this pixel
    }
  }
}

void CleanMapAll(){
  for(auto &it:fImpMap){
    if(it.second->Pin1T()<0) continue; // this pixel is not implanted;
    //TODO::it.second->Fill TTree;
    double pin1e = it.second->Pin1E();
    double tof   = it.second->I2S();
    tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
    h1->Fill(tof, pin1e);
    count1++;
    for(int m=0;m<6;m++){
      if(blob[m]->IsInside(tof,pin1e)){
        h4[m]->Fill(fDecCounter[it.first]);
      }
    }
    it.second->Clear();
    fDecCounter[it.first]=0;
  }
}


void hist(){

  TChannel::ReadDetMapFile();
 
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }
   
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate2.root");
  ban = (TCutG *)cutf2->Get("ban");

  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  cutstrip = (TCutG *)cutf3->Get("cutstrip");
  
  h1 = new TH2D("pid","pid",1000,8000,18000,300,4200,7200);
  for(int i=0;i<6;i++){
    h2[i] = new TH2D(Form("te%i",i)  ,Form("decay time(ms) vs gamma E(keV) in %s",blob[i]->GetName()),1000 ,0,1000,4000,0,4000);
    h3[i] = new TH1D(Form("t%i",i)   ,Form("decay time(ms)in %s",blob[i]->GetName())                 ,10000,0,1000);
    h4[i] = new TH1D(Form("size%i",i),Form("decay size to each imp in %s",blob[i]->GetName())        ,1000 ,0,1000);
  }
  
  for(int runnum = 1040;runnum<1080;runnum++){
    tofpar = TOFCorrection::Get()->ReadFile(runnum);
    tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    if(tofpar.empty()){
      printf("run%i.root doesn't have TOF correction parameters\n",runnum);
      continue;
    }   
 
    TChain *chimp = new TChain("event");
    TChain *chdec = new TChain("event");
    chimp->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/implant/event%i*.root",runnum));
    chdec->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/decay/highgain/event%i*.root",runnum));
    BCSEvent *fimp = new BCSEvent;
    BCSEvent *fdec = new BCSEvent;
    chimp->SetBranchAddress("BCSEvent", &fimp);
    chdec->SetBranchAddress("BCSEvent", &fdec);

    long nimp = chimp->GetEntries();
    long ndec = chdec->GetEntries();
    long ximp = 0;
    long xdec = 0;
    

    InitMap();
    for(ximp=0;ximp<nimp;ximp++){
      fimp->Clear();
      chimp->GetEntry(ximp);  
      bool flag_veto = false;
      if(fimp->LGFSize()==0 || fimp->LGBSize()==0) continue;
      if(fimp->SSSDSize()>0) { flag_veto = veto(fimp->LGFMax(), fimp->SSSD()); }
      if(flag_veto) continue;
      count0++;
      double pin1t = fimp->Pin1T();
      pixel piximp = fimp->LGPixel();
      if(piximp.first<0  || piximp.first>39) continue;
      if(piximp.second<0 || piximp.second>39) continue;
      
      CleanMap(pin1t);     
      //correlation
      while(xdec<ndec){
        fdec->Clear();
        chdec->GetEntry(xdec++);
        if(fdec->HGFMax().GetEnergy()<0 || fdec->HGBMax().GetEnergy()<0) continue;
        double dect = fdec->HGFMax().GetTimestamp();
        if(dect>pin1t) {
          xdec--;
          break;
        }
        pixel pixdec = fdec->HGPixel();
        if(pixdec.first<0  || pixdec.first>39) continue;
        if(pixdec.second<0 || pixdec.second>39) continue;
        if(fImpMap[pixdec]->Pin1T()>0){ 
          double decaytime = dect - fImpMap[pixdec]->Pin1T();
          if(decaytime>1e9 || decaytime<0) continue; // negelect the decay after 1sec and must come after the implant 
          //TODO: correlate this decay to the implant on this pixel;
          fDecCounter[pixdec]++;
          //=== fill gamma hist ===//
          double ctof  = fImpMap[pixdec]->I2S();
          double cpin1e= fImpMap[pixdec]->Pin1E();
          ctof = tofpar[2]*ctof*ctof + tofpar[1]*ctof + tofpar[0]+tofpar1[0];
          decaytime = decaytime/1e6; // unit: ms;
          for(int m=0;m<6;m++){
            if(blob[m]->IsInside(ctof,cpin1e)){
              h3[m]->Fill(decaytime);
            }
          }
          if(fdec->HPGeSize()>0){
            for(int i=0;i<fdec->HPGeSize();i++){
              double gamE = fdec->HPGe()[i].GetEnergy();
              if(gamE<10 || gamE>4000) continue;
              for(int m=0;m<6;m++){
                if(blob[m]->IsInside(ctof,cpin1e)){
                  h2[m]->Fill(decaytime, gamE);
                }
              }
            }
          }
          //=========== gamma info end ==============//
        }
         
      }// correlation end      


      if(fBLMap[piximp]>0){ //reset bl clock;
        fBLMap[piximp] = pin1t;
        count3++;
      }else{
        if(fImpMap[piximp]->Pin1T()>0){ //back to back
          fImpMap[piximp]->Clear();
          fDecCounter[piximp] = 0;
          fBLMap[piximp] = pin1t;
          count2+=2;
        }else{
          fImpMap[piximp]->Clear();
          fImpMap[piximp]->fHits = fimp->fHits;
        }
      }
      
      if((ximp%500)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum, ximp, nimp);
        fflush(stdout);
      }

    }
    
    //to collect decays from the last implant event
    while(xdec<ndec){
      fdec->Clear();
      chdec->GetEntry(xdec++);
      if(fdec->HGFMax().GetEnergy()<0 || fdec->HGBMax().GetEnergy()<0) continue;
      double dect = fdec->HGFMax().GetTimestamp();
      pixel pixdec = fdec->HGPixel();
      if(pixdec.first<0  || pixdec.first>39) continue;
      if(pixdec.second<0 || pixdec.second>39) continue;
      if(fImpMap[pixdec]->Pin1T()>0){ 
        double decaytime = dect - fImpMap[pixdec]->Pin1T();
        if(decaytime>1e9 || decaytime<0) continue; // negelect the decay after 1sec; 
        //TODO: correlate this decay to the implant on this pixel;
        fDecCounter[pixdec]++;
        //=== fill gamma hist ===//
        double ctof  = fImpMap[pixdec]->I2S();
        double cpin1e= fImpMap[pixdec]->Pin1E();
        ctof = tofpar[2]*ctof*ctof + tofpar[1]*ctof + tofpar[0]+tofpar1[0];
        decaytime = decaytime/1e6; // unit: ms;
        for(int m=0;m<6;m++){
          if(blob[m]->IsInside(ctof,cpin1e)){
            h3[m]->Fill(decaytime);
          }
        }
        if(fdec->HPGeSize()>0){
          for(int i=0;i<fdec->HPGeSize();i++){
            double gamE = fdec->HPGe()[i].GetEnergy();
            if(gamE<10 || gamE>4000) continue;
            for(int m=0;m<6;m++){
              if(blob[m]->IsInside(ctof,cpin1e)){
                h2[m]->Fill(decaytime, gamE);
              }
            }
          }
        }
        //=========== gamma info end ==============//
      }
       
    }// correlation end      

    CleanMapAll(); 
    printf("run%i.root: on entry %lu / %lu \n",runnum, ximp, nimp);
  }


   
  std::cout<< "out of veto          = " << count0 << std::endl;
  std::cout<< "imp written in TTree = " << count1 << std::endl;
  std::cout<< "back to back         = " << count2 << std::endl;
  std::cout<< "implanted on blocked = " << count3 << std::endl;
   
  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/correlation3/hist_blobs.root","recreate");
  h1->Write();
  for(int m=0;m<6;m++){
    h2[m]->Write();
    h3[m]->Write();
    h4[m]->Write();
  }
  newf->Close();








}
