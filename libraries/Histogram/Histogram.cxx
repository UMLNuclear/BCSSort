
#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<map>
#include<string>

#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>

#include <Histogram.h>
#include <OutputManager.h>
#include <util.h>
#include <TChannel.h>
#include <TOFCorrection.h>
#include <DetHit.h>
#include <Implant.h>
#include <BCSint.h>
#include <ddaschannel.h>




//======================== Constructor & Destructor========================//


Histogram *Histogram::fHistogram = 0;

Histogram *Histogram::Get(){
  if(!fHistogram){
    fHistogram = new Histogram;
  }
  return fHistogram;
}

Histogram::Histogram(){};
Histogram::~Histogram(){};




//========================== Event File Sort ===========================//
// Get pid and momentum before and after TOF correction//

std::map<int,double[4][4]> Histogram::ReadMat(std::string filename){

  std::map<int, double[4][4]> xtalmat;
  std::ifstream infile;
  std::string line;
  infile.open(filename.c_str());
  int i=0;
  while(getline(infile,line)){
    int det;
    double c0; 
    double c1; 
    double c2; 
    double c3; 

    std::stringstream ss(line);
    ss >> det;
    ss >> c0;
    ss >> c1;
    ss >> c2;
    ss >> c3;

    xtalmat[det][i%4][0] = c0;
    xtalmat[det][i%4][1] = c1;
    xtalmat[det][i%4][2] = c2;
    xtalmat[det][i%4][3] = c3;
    i +=1;
  }

  return xtalmat;
}

void Histogram::EventSort(){
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  std::map<int,double[4][4]> xtalmat = ReadMat();

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long n = gChain->GetEntries();
  //n = 1e5;
  long x = 0;
  int num = 0;
  int numy = 0;
  int numz = 0;
  double adde = 0;

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    int pin1c = 0;
    if(fevent->SSSD().size()>0){
      double sssdEmax = fevent->SSSD().at(0).GetEnergy();
      double sssdEmin = fevent->SSSD().at(0).GetEnergy();
      for(auto &it:fevent->SSSD()){
        if(it.GetEnergy()>sssdEmax) sssdEmax = it.GetEnergy();
        if(it.GetEnergy()<sssdEmin) sssdEmin = it.GetEnergy();
      }
      FillHistogram("eventsize_sssdEmax", 200,0,200,fevent->Size(), 2e3,0,2e4,sssdEmax);
      FillHistogram("eventsize_sssdEmin", 200,0,200,fevent->Size(), 2e3,0,2e4,sssdEmin);
      FillHistogram("Pin1E_sssdEmax",     2e3,0,2e4,fevent->Pin1E(), 2e3,0,2e4,sssdEmax);
      FillHistogram("Pin1E_sssdEmin",     2e3,0,2e4,fevent->Pin1E(), 2e3,0,2e4,sssdEmin);
      FillHistogram("TOF_sssdEmax",     2e3,0,2e4,fevent->I2S(), 2e3,0,2e4,sssdEmax);
      FillHistogram("TOF_sssdEmin",     2e3,0,2e4,fevent->I2S(), 2e3,0,2e4,sssdEmin);

    }
    FillHistogram("event_size",200,0,200,fevent->Size());
    FillHistogram("PID", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E());
    if(fevent->SSSD().size()==0) FillHistogram("PID_nosssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E()); //no SSSD
    else FillHistogram("PID_hassssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E());// has SSSD
    if(fevent->LGFSize()>0 || fevent->LGBSize()>0 || fevent->HGFSize()>0 || fevent->HGBSize()>0){ //hasDSSD
      FillHistogram("event_size_hasDSSD",200,0,200,fevent->Size());
      FillHistogram("PID_hasDSSD", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E());
      //hasDSSD + no SSSD
      if(fevent->SSSD().size()==0) FillHistogram("PID_hasDSSD_nosssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E()); 
      else FillHistogram("PID_hasDSSD_hassssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E()); // hasDSSD + has SSSD
      if(fevent->Pin1E()>0){ //hasDSSD + Pin1>0
        if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
          double dtime = 0;
          if(fevent->LGFSize()>0){
            dtime = fevent->LGFMax().GetTimestamp() - fevent->Pin1T();  
          }else{
            dtime = fevent->LGBMax().GetTimestamp() - fevent->Pin1T();  
          }
          dtime = dtime/1000.;
          FillHistogram("timedif_dssdlopin1", 400,-20,20,dtime, 5e3,0,2e4,fevent->Pin1E());
        }
      }
    }else{//noDSSD
      FillHistogram("event_size_noDSSD",200,0,200,fevent->Size());
      FillHistogram("PID_noDSSD", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E());
      //noDSSD + no SSSD
      if(fevent->SSSD().size()==0) FillHistogram("PID_noDSSD_nosssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E());
      else FillHistogram("PID_noDSSD_hassssd", 2e3,0,2e4,fevent->I2S(),5e3,0,20e3,fevent->Pin1E()); //noDSSD + has SSSD
      if(fevent->HPGe().size()==0 && fevent->LaBr().size()==0){ //no DSSD + not Clover Only
        for(auto &it:fevent->fHits){
          FillHistogram("summary_mystery",2e3,0,2e4,it.GetEnergy(), 300,0,300,it.GetNumber());
          FillHistogram("sumsize_mystery",200,0,200,fevent->Size(), 300,0,300,it.GetNumber());
        }
      }
    }
    for(auto &it:fevent->fHits){
      FillHistogram("summary", 2000,0,20000,it.GetEnergy(),300,0,300,it.GetNumber());
      FillHistogram("sum_size", 200,0,200,fevent->Size(), 300,0,300,it.GetNumber());
      if(fevent->LGFSize()>0 || fevent->LGBSize()>0 || fevent->HGFSize()>0 || fevent->HGBSize()>0){
        if(fevent->Pin1E()>0 && (fevent->LGFSize()>0 || fevent->LGBSize()>0) && it.GetNumber()==181) pin1c++;
        FillHistogram("summary_hasDSSD", 2000,0,20000,it.GetEnergy(),300,0,300,it.GetNumber());
        FillHistogram("sum_size_hasDSSD", 200,0,200,fevent->Size(), 300,0,300,it.GetNumber());
      }else{
        FillHistogram("summary_noDSSD", 2000,0,20000,it.GetEnergy(),300,0,300,it.GetNumber());
        FillHistogram("sum_size_noDSSD", 200,0,200,fevent->Size(), 300,0,300,it.GetNumber());
      }
      for(auto &it1:fevent->fHits){
        double dt = it.GetTimestamp() - it1.GetTimestamp();
        dt = dt/1000.;
        FillHistogram("timedif_hits",400,-20,20,dt);
        if(fevent->LGFSize()>0 || fevent->LGBSize()>0 || fevent->HGFSize()>0 || fevent->HGBSize()>0){
          FillHistogram("timedif_hits_hasDSSD",400,-20,20,dt);
          if(fevent->Pin1E()>0){ //Implant
            FillHistogram("timedif_Implant",400,-20,20,dt);
            if(it.GetNumber()>180 && it.GetNumber()<184 && it1.GetNumber()>180 && it1.GetNumber()<184){
              FillHistogram("timedif_Implant_pin",400,-20,20,dt);
            }
            if((it.GetNumber()>39 && it.GetNumber()<80) || (it.GetNumber()>119 && it.GetNumber()<160)){
              if((it1.GetNumber()>39 && it1.GetNumber()<80) || (it1.GetNumber()>119 && it1.GetNumber()<160)){
                FillHistogram("timedif_Implant_dssdlo",400,-20,20,dt);
              }
            }
            if((it.GetNumber()>39 && it.GetNumber()<80) || (it.GetNumber()>119 && it.GetNumber()<160) 
            || (it.GetNumber()>180 && it.GetNumber()<184)){
              if((it1.GetNumber()>39 && it1.GetNumber()<80) || (it1.GetNumber()>119 && it1.GetNumber()<160) 
              || (it1.GetNumber()>180 && it1.GetNumber()<184)){
                FillHistogram("timedif_Implant_pindssdlo",400,-20,20,dt);
              }
            }
            if(it.GetNumber()<160 && it1.GetNumber()<160){
              FillHistogram("timedif_Implant_dssd",400,-20,20,dt);
            }
            if(it.GetNumber()>159 && it.GetNumber()<176 && it1.GetNumber()>159 && it1.GetNumber()<176){
              FillHistogram("timedif_Implant_sssd",400,-20,20,dt);
            }
          }else{//Decay
            FillHistogram("timedif_Decay",400,-20,20,dt);
            if(it.GetNumber()<160 && it1.GetNumber()<160){
              FillHistogram("timedif_Decay_dssd",400,-20,20,dt);
            }
            if(it.GetNumber()>159 && it.GetNumber()<176 && it1.GetNumber()>159 && it1.GetNumber()<176){
              FillHistogram("timedif_Decay_sssd",400,-20,20,dt);
            }
          }
        }else{
          FillHistogram("timedif_hits_noDSSD",400,-20,20,dt);
        }
      }
    }
    FillHistogram("pin1size_pin1lo",50,0,50,pin1c);


    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("event_outputdo%s.root",runnum.c_str()));
  //SaveHistograms("eventdo0048-00.root");
  return;


}

void Histogram::EventHitPad(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());

  std::vector<TCutG* >pidcuts;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/pidcut_48.root");
  TIter keys1 (mycut1->GetListOfKeys());  
  while(TKey *key1 = (TKey*)keys1.Next()){
    pidcuts.push_back((TCutG*)key1->ReadObj());
  }

  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  long nentries = gChain->GetEntries();
  //nentries = 1e1;
  long x = 0;
  double dt = 0;

  for(x=0;x<nentries;x++){
    gChain->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
        dt = fevent->LGFMax().GetTimestamp() - fevent->LGBMax().GetTimestamp();
        dt = dt/1000.;
        FillHistogram("imp_dt_loFandB", 2000,-10,10,dt);
      }
      if(fevent->DSSDloT()>0){
        dt = fevent->DSSDloT() - fevent->Pin1T();
        dt = dt/1000.;
        FillHistogram("imp_pin1E_dt", 2000,-10,10,dt, 2e3,0,2e4,fevent->Pin1E());
        FillHistogram("imp_TOF_dt", 2000,-10,10,dt, 2e3,0,2e4,fevent->I2S());
      }
    }else{ // no Pin1
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        dt = fevent->HGFMax().GetTimestamp() - fevent->HGBMax().GetTimestamp();
        dt = dt/1000.;
        FillHistogram("dec_dt_hiFandB", 2000,-10,10,dt);
      }
      if(fevent->DSSDhiT()>0){
        for(auto &it:fevent->HPGe()){
          dt = fevent->DSSDhiT() - it.GetTimestamp();
          dt = dt/1000.;
          FillHistogram("dec_HPGeE_dt", 2000,-10,10,dt, 8e3,10,4e3,it.GetEnergy());
        }
        for(auto &it:fevent->LaBr()){
          dt = fevent->DSSDhiT() - it.GetTimestamp();
          dt = dt/1000.;
          FillHistogram("dec_LaBrE_dt", 2000,-10,10,dt, 8e3,10,4e3,it.GetEnergy());
        }
      }
    }

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,nentries);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,nentries);
  SaveHistograms(Form("timedif_eventdo%s.root",runnum.c_str()));
  return;
}




//=========================== ListSort ===================================//

void Histogram::Process(std::vector<DetHit> vec, std::map<int,double[4][4]> xtalmat){
  int sssd_flag = 0;
  std::vector<DetHit> FL;
  std::vector<DetHit> BL;
  std::vector<DetHit> Ge;
  DetHit pin1_i2s;
  DetHit pin1;
  DetHit i2s_i2n;
  int num = 0;
  int numy = 0;
  int numz = 0;
  for(size_t y=0;y<vec.size();y++){
    FillHistogram("Summary",16e3,0,16e3,vec[y].GetCharge(), 300,0,300,vec[y].GetNumber());
    FillHistogram("Summary_cal",16e3,0,16e3,vec[y].GetEnergy(), 300,0,300,vec[y].GetNumber());
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
        break;

      case 208 ... 271: //HPGe
        Ge.push_back(vec[y]);
        break;

      default:
        break;
    }
  }
  if(pin1.GetCharge()>100){
    FillHistogram("PID",2e3,0,2e4,pin1_i2s.GetCharge(),10e3,0,20e3,pin1.GetCharge());
    FillHistogram("Pin1E",10e3,0,20e3,pin1.GetCharge());
  }
  std::map<int,Clover> clmap; 
  if(Ge.size()>0){
    for(auto it:Ge){
      num = it.GetNumber()-208;
      if(it.GetEnergy()>10 && it.GetEnergy()<4000){
        int clnum = num/4;
        int xtalnum = num%4;
        if(clnum == 8) {FillHistogram(Form("single_cl8_%i",xtalnum),8e3,0,4e3,it.GetEnergy());}
        if(clnum == 9) {FillHistogram(Form("single_cl9_%i",xtalnum),8e3,0,4e3,it.GetEnergy());} 
        //if(it.GetEnergy()>1450)
        //printf("channum = [%i][%i] \t Charge = %f \t Energy = %f\n\n",clnum,xtalnum,it.GetCharge(),it.GetEnergy());
        clmap[clnum].Add(it);
        FillHistogram("single",8000,0,4000,it.GetEnergy());
        FillHistogram("single_char",8000,0,4000,it.GetCharge());
        FillHistogram(Form("single_%i",clnum),8000,0,4000,it.GetEnergy());
      }  
    }

    std::map<int,Clover>::iterator it;
    for(it=clmap.begin();it!=clmap.end();it++){
      it->second.SetAddE();
      double adde = 0;
      for(size_t m=0;m<it->second.Size();m++){
        adde += it->second.fXtal[m].GetEnergy();
        for(size_t z=0;z<it->second.Size();z++){
          numy = (it->second.fXtal[m].GetNumber()-208) % 4;
          numz = (it->second.fXtal[z].GetNumber()-208) % 4;
          adde += xtalmat[it->first][numy][numz]*it->second.fXtal[z].GetEnergy();
        }
      }
      FillHistogram("Addback",8000,0,4000,it->second.AddbackE());
      FillHistogram(Form("Addback_%i",it->first),8000,0,4000,it->second.AddbackE());
      it->second.SetAddE(adde);
      FillHistogram("Addback_ct",8000,0,4000,it->second.AddbackE());
      FillHistogram(Form("Addback_ct_%i",it->first),8000,0,4000,it->second.AddbackE()); 
    }
  }
}


void Histogram::ListSort(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  int faddress;
  int fnumber;
  double ftimestamp;
  double fcharge;
  gChain->SetBranchAddress("address",   &faddress);
  gChain->SetBranchAddress("number",    &fnumber);
  gChain->SetBranchAddress("timestamp", &ftimestamp);
  gChain->SetBranchAddress("charge",    &fcharge);

  TChannel::ReadDetMapFile();
  std::map<int,double[4][4]> xtalmat = ReadMat();

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
    if((fhit->GetTimestamp()-starttime) > buildtime){
      starttime = fhit->GetTimestamp();
      Process(vec_hit, xtalmat);
      vec_hit.clear();
    }
    vec_hit.push_back(*fhit);
    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("list_output%s.root",runnum.c_str()));
  return;

}



//================================ LogX  ================================//
void Histogram::LogX(){
  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta",&fbeta);
  TChannel::ReadDetMapFile();

  int Nbins = 1000;
  double xlow = 4; // unit: nanosecond
  double xhigh = 1e9; //5e8 ns = 500ms
  double dx = log(xhigh/xlow)/(double)Nbins;
  Double_t edges[Nbins+1];
  edges[0] = 0;
  for(int i=0;i<Nbins+1;i++){
    edges[i+1] = exp(log(xlow) + (double)i*dx);
  }

  long x=0;
  long n = gChain->GetEntries();

  std::vector<double> tofpara;
  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;
  if(tofpara.empty()) return;

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newsubblobNa.root"); // subcuts from pid only Na32 after TOF fluctuation.
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  std::vector<TCutG *>cut2;
  TFile *mycut2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/decaytimecuts.root"); // cuts from decaytime and HPGe.Energy; timecuts[0] is "LHE" cut
  TIter keys2 (mycut2->GetListOfKeys());
  while(TKey *key2 = (TKey*)keys2.Next()) {
    cut2.push_back((TCutG*)key2->ReadObj());
  }
  TH1D *hist = new TH1D("Logx_dtime", "LogX_dtime", Nbins,edges);
  TH1D *hist1 = new TH1D("Logx_dtime1", "LogX_dtime1", Nbins,edges);
  TH1D *hist2 = new TH1D("Logx_dtime2", "LogX_dtime2", Nbins,edges);
  TH1D *hist3 = new TH1D("Logx_dtime3", "LogX_dtime3", Nbins,edges);
  TH1D *hist4 = new TH1D("Logx_dtime4", "LogX_dtime4", Nbins,edges);

  TH1D *ht = new TH1D("Logx_decaytime", "LogX_decaytime", Nbins,edges);
  TH1D *ht1 = new TH1D("Logx_decaytime1", "LogX_decaytime1", Nbins,edges);
  TH1D *ht2 = new TH1D("Logx_decaytime2", "LogX_decaytime2", Nbins,edges);
  TH1D *ht3 = new TH1D("Logx_decaytime3", "LogX_decaytime3", Nbins,edges);
  TH1D *ht4 = new TH1D("Logx_decaytime4", "LogX_decaytime4", Nbins,edges);

  TFile *newf = new TFile(Form("logdecay_%s.root",runnum.c_str()),"recreate");

  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    double charge = fbeta->fImplant.fI2S;
    charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
    if(fbeta->fImplant.Stopped()){
      int count = 0;
      for(auto &it1: fbeta->fDecay){
        count++;
        for(auto &it2: it1.fGe){
          double dt = it1.GetTimestamp() - it2.GetTimestamp();
          for(size_t m=0;m<cut1.size();m++){
            if(m==1) {
              if(cut1[m]->IsInside(charge ,fbeta->fImplant.fPIN1E)){
                double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
                dtime = fabs(dtime); //unit: ns;
                if(count<2){
                  hist1->Fill(dtime);
                }
                if(count<3){
                  hist2->Fill(dtime);
                }
                if(count<4){
                  hist3->Fill(dtime);
                }
                if(count<5){
                  hist4->Fill(dtime);
                }
                hist->Fill(dtime);

                if(cut2[0]->IsInside(dt, it2.GetEnergy())){
                  if(count<2){
                    ht1->Fill(dtime);
                  }
                  if(count<3){
                    ht2->Fill(dtime);
                  }
                  if(count<4){
                    ht3->Fill(dtime);
                  }
                  if(count<5){
                    ht4->Fill(dtime);
                  }
                  ht->Fill(dtime);
                }


              }
            }          
          }

        }
      }
    }
    if((x%20000)==0){
      printf("on entry %lu / %lu  \r",x,n);
      fflush(stdout);
    }
  }
  printf("      on entry %lu   /   %lu    \n",x,n);
  hist1->Write();
  hist2->Write();
  hist3->Write();
  hist4->Write();
  hist->Write();
  ht->Write();
  ht1->Write();
  ht2->Write();
  ht3->Write();
  ht4->Write();
  newf->Close();
  return;
}


//========================================== ================================//






