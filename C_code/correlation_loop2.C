
#include<Implant.h>
#include<TOFCorrection.h>
#include<globals.h>
#include<util.h>


// global variables;
TChain *chimp = new TChain("implant");
TChain *chdec = new TChain("decay");

Implant *fimp = new Implant;
Decay *fdec = new Decay;


std::map<pixel,Implant*> fImpMap;
std::map<pixel,std::vector<Decay> *> fDecMap;
std::map<pixel,std::vector<Decay> *> fDecMap2;
std::map<pixel,std::vector<Decay> *> fDecMap3;
std::map<pixel,std::vector<Decay> *> fDecMap4;
std::map<pixel,std::vector<Decay> *> fDecMap5;
std::map<pixel,std::vector<Decay> *> fDecMap6;
std::map<pixel, double> fImpTime;

//==========Initial Map=========//
void Initial(){
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new Implant;
      fDecMap[pix] = new std::vector<Decay>;
      fDecMap2[pix] = new std::vector<Decay>;
      fDecMap3[pix] = new std::vector<Decay>;
      fDecMap4[pix] = new std::vector<Decay>;
      fDecMap5[pix] = new std::vector<Decay>;
      fDecMap6[pix] = new std::vector<Decay>;
      fImpTime[pix] = -1;
    }
  }
}


//====================== Main() =====================//

void correlation_loop(){
  int runnum = 1032;
  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pidcut104.root"); // PID
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }
  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/prompt_Imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); 
  veccut.push_back(pin1E_dt);   

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");

  TChannel::ReadDetMapFile();

  long nimp, ndec;
  long ximp = 0;
  long xdec = 0;
  int c  = 0;
  int c1 = 0; 
  int c2 = 0; 
  int c3 = 0; 
  int c4 = 0; 
  int c5 = 0; 
  int c6 = 0; 
  int c7 = 0; 
  int c8 = 0; 

  double dt;

  chimp->Add(Form("data/5us_tofcor/implant/implant%i*",runnum));
  chdec->Add(Form("data/5us_tofcor/decay/dec_prompt/decay%i*",runnum));
  std::vector<double> tofpara = TOFCorrection::Get()->ReadFile(runnum, "/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
  if(tofpara.empty()){
    printf("run%i.root file doesn't have tof parameters\n", runnum);
    return;
  }
  chimp->SetBranchAddress("Implant", &fimp);
  chdec->SetBranchAddress("Decay", &fdec);
  nimp = chimp->GetEntries();
  ndec = chdec->GetEntries();

  Initial();

  double oldT;
  double twindow = 50;
  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    fimp->fI2S = fimp->fI2S + tofpara[0];
    pixel piximp = fimp->GetPixel();
    if(piximp.first<0 ||  piximp.second<0) continue;
    if(piximp.first>39 ||  piximp.second>39) continue;
    if(fImpTime[piximp]>0){
      while(xdec<ndec){
        chdec->GetEntry(xdec++);
        pixel pixdec = fdec->GetPixel();
        if(fdec->DSSDhiT()>fimp->DSSDloT()){
          xdec--;
          Implant *current = fImpMap[piximp];
          std::vector<Decay> *vdec = fDecMap[piximp];
          std::vector<Decay> *vdec2 = fDecMap2[piximp];
          std::vector<Decay> *vdec3 = fDecMap3[piximp];
          std::vector<Decay> *vdec4 = fDecMap4[piximp];
          std::vector<Decay> *vdec5 = fDecMap5[piximp];
          std::vector<Decay> *vdec6 = fDecMap6[piximp];
          dt = current->DSSDloT() - current->fPIN1T;
          dt = dt/1000.; //us
          if(veccut[0]->IsInside(dt, current->fI2S) && veccut[1]->IsInside(dt,current->fPIN1E)){
            double i2 = current->fI2S_I2N;
            double tof = current->fI2S;
            if((Ne->IsInside(tof, current->fPIN1E)) && (i2>10800 && i2<13800) && (tof>11100 && tof<13000)){
              for(auto &it1:*vdec){
                if(it1.fDecayTime>50) continue;
                for(auto &it11:it1.fGe){
                  if(it11.GetEnergy()<10 || it11.GetEnergy()>4000) continue;
                  FillHistogram("dt_singles_1pixel",500,0,50,it1.fDecayTime, 4000,0,4000,it11.GetEnergy());
                } 
              }
              for(auto &it2:*vdec2){
                if(it2.fDecayTime>50) continue;
                for(auto &it22:it2.fGe){
                  if(it22.GetEnergy()<10 || it22.GetEnergy()>4000) continue;
                  FillHistogram("dt_singles_9pixel",500,0,50,it2.fDecayTime, 4000,0,4000,it22.GetEnergy());
                } 
              }
              for(auto &it3:*vdec3){
                if(it3.fDecayTime>50) continue;
                for(auto &it33:it3.fGe){
                  if(it33.GetEnergy()<10 || it33.GetEnergy()>4000) continue;
                  FillHistogram("dt_singles_25pixel",500,0,50,it3.fDecayTime, 4000,0,4000,it33.GetEnergy());
                } 
              }
              for(auto &it4:*vdec4){
                if(it4.fDecayTime>50) continue;
                for(auto &it44:it4.fGe){
                  if(it44.GetEnergy()<10 || it44.GetEnergy()>4000) continue;
                  FillHistogram("dt_singles_49pixel",500,0,50,it4.fDecayTime, 4000,0,4000,it44.GetEnergy());
                } 
              }
              for(auto &it5:*vdec5){
                if(it5.fDecayTime>50) continue;
                for(auto &it55:it5.fGe)l.first+ext;i++){
            if(i<0 || i>39) continue;
            for(int j=pixdec.second-ext;j<=pixdec.second+ext;j++){
              if(j<0 || j>39) continue;
              pixel temp = std::make_pair(i,j);
              if(fImpTime[temp]>0){
                if(fdec->DSSDhiT()<fImpTime[temp]) continue;
                double tdif = fdec->DSSDhiT() - fImpMap[temp]->fPIN1T;
                if(tshort>tdif){
                  tshort = tdif;
                  best = temp;
                }
              }
            }
          }
          if(best.first<0 || best.second<0) continue;
          dt = fdec->DSSDhiT() - fImpMap[best]->fPIN1T;
          dt = dt/1000000.;//unit:ms 
          fdec->SetDecayTime(dt);
          if(ext==0) fDecMap[best]->push_back(*fdec);
          if(ext==1) fDecMap2[best]->push_back(*fdec);
          if(ext==2) fDecMap3[best]->push_back(*fdec);
          if(ext==3) fDecMap4[best]->push_back(*fdec);
          if(ext==4) fDecMap5[best]->push_back(*fdec);
          if(ext==5) fDecMap6[best]->push_back(*fdec);
        }
      }
    }
    fimp->Copy(*fImpMap[piximp]);
    fImpTime[piximp] = fimp->DSSDloT();

    if((ximp%5000)==0){
      printf("  on entry %lu / %lu    \r", ximp, nimp);
      fflush(stdout);
    }
  }


  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel piximp = std::make_pair(m,n);
      if(fImpTime[piximp]>0){
        std::vector<Decay> *vdec = fDecMap[piximp];
        std::vector<Decay> *vdec2 = fDecMap2[piximp];
        std::vector<Decay> *vdec3 = fDecMap3[piximp];
        std::vector<Decay> *vdec4 = fDecMap4[piximp];
        std::vector<Decay> *vdec5 = fDecMap5[piximp];
        std::vector<Decay> *vdec6 = fDecMap6[piximp];
        dt = fImpMap[piximp]->DSSDloT() - fImpMap[piximp]->fPIN1T;
        dt = dt/1000.; //us
        if(veccut[0]->IsInside(dt, fImpMap[piximp]->fI2S) && veccut[1]->IsInside(dt,fImpMap[piximp]->fPIN1E)){
          double i2 = fImpMap[piximp]->fI2S_I2N;
          double tof = fImpMap[piximp]->fI2S;
          if((Ne->IsInside(tof, fImpMap[piximp]->fPIN1E)) && (i2>10800 && i2<13800) && (tof>11100 && tof<13000)){
            for(auto &it1:*vdec){
              if(it1.fDecayTime>50) continue;
              for(auto &it11:it1.fGe){
                if(it11.GetEnergy()<10 || it11.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_1pixel",500,0,50,it1.fDecayTime, 4000,0,4000,it11.GetEnergy());
              } 
            }
            for(auto &it2:*vdec2){
              if(it2.fDecayTime>50) continue;
              for(auto &it22:it2.fGe){
                if(it22.GetEnergy()<10 || it22.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_9pixel",500,0,50,it2.fDecayTime, 4000,0,4000,it22.GetEnergy());
              } 
            }
            for(auto &it3:*vdec3){
              if(it3.fDecayTime>50) continue;
              for(auto &it33:it3.fGe){
                if(it33.GetEnergy()<10 || it33.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_25pixel",500,0,50,it3.fDecayTime, 4000,0,4000,it33.GetEnergy());
              } 
            }
            for(auto &it4:*vdec4){
              if(it4.fDecayTime>50) continue;
              for(auto &it44:it4.fGe){
                if(it44.GetEnergy()<10 || it44.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_49pixel",500,0,50,it4.fDecayTime, 4000,0,4000,it44.GetEnergy());
              } 
            }
            for(auto &it5:*vdec5){
              if(it5.fDecayTime>50) continue;
              for(auto &it55:it5.fGe){
                if(it55.GetEnergy()<10 || it55.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_81pixel",500,0,50,it5.fDecayTime, 4000,0,4000,it55.GetEnergy());
              } 
            }
            for(auto &it6:*vdec6){
              if(it6.fDecayTime>50) continue;
              for(auto &it66:it6.fGe){
                if(it66.GetEnergy()<10 || it66.GetEnergy()>4000) continue;
                FillHistogram("dt_singles_121pixel",500,0,50,it6.fDecayTime, 4000,0,4000,it66.GetEnergy());
              } 
            }
          }
        }
      }
    }
  }



  SaveHistograms(Form("correlation_loop2%i.root",runnum));
  //std::cout << "imp hit again within time window = " << c1 << std::endl;


  return;
}
