
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
  TCutG *Na = (TCutG *)cutf1->Get("Na");

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

  double twindow = 100;
  for(ximp=0;ximp<nimp;ximp++){
    chimp->GetEntry(ximp);
    fimp->fI2S = fimp->fI2S + tofpara[0];
    pixel piximp = fimp->GetPixel();
    if(piximp.first<0 ||  piximp.second<0) continue;
    if(piximp.first>39 ||  piximp.second>39) continue;
    dt = fimp->DSSDloT() - fimp->fPIN1T;
    dt = dt/1000.; //us
    if(veccut[0]->IsInside(dt, fimp->fI2S) && veccut[1]->IsInside(dt,fimp->fPIN1E)){
      double i2 = fimp->fI2S_I2N;
      double tof = fimp->fI2S;
      if((Ne->IsInside(tof, fimp->fPIN1E)) && (i2>10800 && i2<13800) && (tof>11100 && tof<13000)){
        while(xdec<ndec){
          chdec->GetEntry(xdec++);
          if(fdec->DSSDhiT()<fimp->fPIN1T) continue;
          if(fdec->DSSDhiT()>(fimp->fPIN1T+twindow*1e6)){
            xdec--;
            break;
          }
          for(auto &it:fdec->fGe){
            if(it.GetEnergy()<10 || it.GetEnergy()>4000) continue;
            double tdif = fdec->DSSDhiT() - fimp->fPIN1T;
            tdif = tdif/1000000.;
            FillHistogram("dt_singles_nopixel",500,0,50,tdif,4000,0,4000,it.GetEnergy());
          }
        }
      }    
    }


      if((ximp%5000)==0){
        printf("  on entry %lu / %lu    \r", ximp, nimp);
        fflush(stdout);
      }
    }
    printf("  on entry %lu / %lu    \n", ximp, nimp);
    SaveHistograms(Form("correlation_loop1%i",runnum));

    //std::cout<<"total imp = "      <<nimp<<std::endl;
    //std::cout<<"time window(ms) = "      <<twindow<<std::endl; 
    //std::cout<<"dt imp shorter than time window = "  <<c1<<std::endl;
    //std::cout<<"the former one is 31Ne = "           <<c3<<std::endl;
    //std::cout<<"the latter one is 31Ne = "           <<c2<<std::endl;



    return;
  }
