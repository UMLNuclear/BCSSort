
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
std::map<pixel, double> fImpTime;

//==========Initial Map=========//
void Initial(){
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      fImpMap[pix] = new Implant;
      fDecMap[pix] = new std::vector<Decay>;
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
  int c1 = 0; // how many implants around one decay
  int c2 = 0; // how many valid implants around one decay;
  int c3 = 0; // the correlated implant has different pixel;
  int c4 = 0; // the correlated valid implant has different pixel;
  int c5 = 0; // how many good implants;
  int c6 = 0; // how many good valid implants;
  int c7 = 0; // how many good valid implants has decay within 50ms;
  int c8 = 0; // how many good valid implants has decay with different pixels within 50ms;

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

  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    fimp->fI2S = fimp->fI2S + tofpara[0];
    pixel piximp = fimp->GetPixel();
    if(piximp.first<0 ||  piximp.second<0) continue;
    if(piximp.first>39 ||  piximp.second>39) continue;
    if(fImpTime[piximp]>0){ //that pixel is implanted;
      while(xdec<ndec){
        chdec->GetEntry(xdec++);
        pixel pixdec = fdec->GetPixel();
        if(pixdec.first<0 || pixdec.second<0) continue;
        if(fdec->DSSDhiT() > fimp->DSSDloT()) {
          bool flag1 = true;
          bool flag = true;
          c3 = 0;
          c4 = 0;
          xdec--;
          Implant *current = fImpMap[piximp];
          std::vector<Decay> *vdec = fDecMap[piximp];
          dt = current->DSSDloT() - current->fPIN1T;
          dt = dt/1000.; //us
          if(veccut[0]->IsInside(dt, current->fI2S) && veccut[1]->IsInside(dt,current->fPIN1E)){
            double i2 = current->fI2S_I2N;
            c5++;
            if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
              c6++;
            }
            for(auto &it:*vdec){
              if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
                if(it.fDecayTime<50){
                  c3++;
                  if(flag){
                    c7++;
                    flag = false;
                  }
                }
              }
              if(it.GetPixel().first!=current->GetPixel().first || it.GetPixel().second!=current->GetPixel().second){
                if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
                  if(it.fDecayTime<50){
                    c4++;
                    if(flag1){
                      c8++;
                      flag1 = false;
                    }
                  }
                }
              }
            }
          }
          current->Clear();
          vdec->clear();
          FillHistogram("Ne31_decaysize",100,0,100,c3);
          FillHistogram("Ne31_DIFdecaysize",100,0,100,c4);
          break;
        }
        double tshort = DBL_MAX;
        pixel best = std::make_pair(-1,-1);
        c1 = 0;// check how many implants around one decay;
        c2 = 0;// check how many valid implants around one decay;(valid = 31Ne + 50ms decaytime)
        for(int i=pixdec.first-1;i<=pixdec.first+1;i++){
          if(i<0 || i>39) continue;
          for(int j=pixdec.second-1;j<=pixdec.second+1;j++){
            if(j<0 || j>39) continue;
            pixel temp = std::make_pair(i,j);
            if(fImpTime[temp]>0 && fImpTime[temp]<fdec->DSSDhiT()){
              c1++;
              double tdif = fdec->DSSDhiT() - fImpMap[temp]->fPIN1T;
              tdif = tdif/1000000.;
              if(tdif<50){
                Implant *tempimp = fImpMap[temp];
                double i2 = tempimp->fI2S_I2N;
                if((Ne->IsInside(tempimp->fI2S, tempimp->fPIN1E)) && (i2>10800 && i2<13800)){
                  c2++;
                }
              }
              if(tshort>fabs(fdec->DSSDhiT() - fImpTime[temp])){
                tshort = fabs(fdec->DSSDhiT() - fImpTime[temp]);
                best = temp;
              }
            }
          }
        }
        if(best.first<0 || best.first>39) continue;
        if(best.second<0 || best.second>39) continue;
        dt = fdec->DSSDhiT() - fImpMap[pixdec]->fPIN1T;
        dt = dt/1000000.;//unit:ms 
        fdec->SetDecayTime(dt);
        fDecMap[best]->push_back(*fdec);
        FillHistogram("ImpSizeofADec_50ms",9,0,9,c1);
        FillHistogram("ValidImpSizeofADec_50ms",9,0,9,c2);
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
      pixel pix = std::make_pair(m,n);
      if(fImpTime[pix]<0) continue;
      Implant *current = fImpMap[pix];
      std::vector<Decay> *vdec = fDecMap[pix];
      dt = current->DSSDloT() - current->fPIN1T;
      dt = dt/1000.; //us
      bool flag1 = true;
      bool flag = true;
      c3 = 0;
      c4 = 0;
      if(veccut[0]->IsInside(dt, current->fI2S) && veccut[1]->IsInside(dt,current->fPIN1E)){
        double i2 = current->fI2S_I2N;
        c5++;
        if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
          c6++;
        }
        for(auto &it:*vdec){
          if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
            if(it.fDecayTime<50){
              c3++;
              if(flag){
                c7++;
                flag = false;
              }
            }
          }
          if(it.GetPixel().first!=current->GetPixel().first || it.GetPixel().second!=current->GetPixel().second){
            if(Ne->IsInside(current->fI2S, current->fPIN1E) && (i2>10800 && i2<13800)){
              if(it.fDecayTime<50){
                c4++;
                if(flag1){
                  c8++;
                  flag1 = false;
                }
              }
            }
          }
        }
      }
      current->Clear();
      vdec->clear();
      FillHistogram("Ne31_decaysize",100,0,100,c3);
      FillHistogram("Ne31_DIFdecaysize",100,0,100,c4);
    }
  }

  SaveHistograms(Form("correlation_loop%i.root",runnum));
  std::cout<<"imp total entries = "          << nimp <<std::endl;
  std::cout<<"dec total entries = "          << ndec <<std::endl;
  std::cout<<"good implant = "               << c5   <<std::endl;
  std::cout<<"good 31Ne implant = "          << c6   <<std::endl;
  //std::cout<<"decay correlated to good imp has different pixels = " << c3   <<std::endl;
  //std::cout<<"decay correlayed to good 31Ne has different pixels = " << c4   <<std::endl;
  std::cout<<"31Ne has decays within 50ms = " << c7   <<std::endl;
  std::cout<<"31Ne has decay with different pixels within 50ms = " << c8   <<std::endl;



  return;
}
