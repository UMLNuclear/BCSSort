

#include <cstdio>
#include <string>

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
void Histogram::EventSort(){


  std::vector<double> tofpara;
  BCSEvent *fevent = new BCSEvent;
  gChain->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  //OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;
  if(tofpara.empty()) return;

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/pidcuts_momen.root");
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  //====================== parameters about TOF and momentum ======================//

  std::vector<double> fitpary1;
  fitpary1.push_back(0.2); 
  fitpary1.push_back(0.22); 
  fitpary1.push_back(0.26); 
  fitpary1.push_back(0.379711); 
  fitpary1.push_back(0.241714); 
  fitpary1.push_back(0.183131); 
  fitpary1.push_back(0.342038); 

  std::vector<double> fitpary2;
  for(size_t i=0;i<fitpary1.size();i++){
    fitpary2.push_back(-fitpary1[i]);
  }

  //==============================================================================//


  long n = gChain->GetEntries();
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(fevent->Pin1E()>0 && fevent->I2S()>0){
      double charge = fevent->I2S();
      charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
      FillHistogram("pid", 2e3,0,2e4,fevent->I2S(), 4e3,0,8e3,fevent->Pin1E());
      FillHistogram("pid_after", 2e3,0,2e4,charge, 4e3,0,8e3,fevent->Pin1E());
      if(fevent->SSSDSize()==0){
        FillHistogram("pid_nosssd", 2e3,0,2e4,fevent->I2S(), 4e3,0,8e3,fevent->Pin1E());
        FillHistogram("pid__nosssd_after", 2e3,0,2e4,charge, 4e3,0,8e3,fevent->Pin1E());
      }
    }
    if(fevent->I2N_I2S()>0 && fevent->Pin1E()>0 && fevent->LGBSize() && fevent->LGFSize()){
    }

    if((x%200000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("event_output%s.root",runnum.c_str())); 

  return;
}



//======================== Beta File Sort =========================//
void  Histogram::BetaSort(){
  std::vector<double> tofpara;
  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta", &fbeta);
  TChannel::ReadDetMapFile();

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);
  tofpara = TOFCorrection::Get()->ReadFile(std::stoi(runnum));
  //OutputManager::Get()->Set(runnum);
  std::cout << "created output tree file: " << runnum << std::endl;
  if(tofpara.empty()) return;

  std::vector<TCutG *> cut1;
  //TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newsubblobNa_024.root"); // subcuts from pid only Na32 after TOF fluctuation and TOF&I2 correction.
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/newsubblobNe.root"); // subcuts from pid only Na32 after TOF fluctuation.
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

  //====================== parameters about TOF and momentum ======================//

  //std::vector<double> fitpary1;
  //fitpary1.push_back(-0.0243118); 
  ////fitpary1.push_back(0.0142666); 
  ////fitpary1.push_back(0.00555406); 
  ////fitpary1.push_back(-0.132843); 

  //std::vector<double> fitpary2;
  //for(size_t i=0;i<fitpary1.size();i++){
  //  fitpary2.push_back(-fitpary1[i]);
  //}




  long n = gChain->GetEntries();
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    double charge = fbeta->fImplant.fI2S;
    charge = tofpara[0]+tofpara[1]*charge+tofpara[2]*charge*charge;
    if(fbeta->fImplant.fPIN1E>0 && fbeta->fImplant.fI2S>0){
      FillHistogram("pid",2e3,0,2e4,fbeta->fImplant.fI2S, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
      FillHistogram("pid_after",2e3,0,2e4,charge, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
    }
    if(fbeta->fImplant.Stopped()){
      bool flag1=true;
      bool flag2=true;

      //double new_y2 = charge +  fitpary2[0]*fbeta->fImplant.fI2S_I2N;        
      //double new_y2_1 = charge +  fitpary2[1]*fbeta->fImplant.fI2S_I2N;        
      if(fbeta->fImplant.fPIN1E>0 && fbeta->fImplant.fI2S>0){
        FillHistogram("pid_nosssd_after",2e3,0,2e4,charge, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
        FillHistogram("pid_nosssd",2e3,0,2e4,fbeta->fImplant.fI2S, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
      }
      //for(size_t m=0;m<fitpary1.size();m++){
      //double new_y = charge +  fitpary2[m]*fbeta->fImplant.fI2S_I2N;        
      //FillHistogram(Form("pid_y2_%.3f",fabs(fitpary1[m])),3e3,0,3e4, new_y, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
      if(fbeta->DecaySize()>0){
        FillHistogram("pid_De",3e3,0,3e4,charge, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
        for(auto &it1:fbeta->fDecay){
          if(it1.fImplantTime<100){
            if(flag1){
              FillHistogram("pid_De100",3e3,0,3e4,charge, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
              flag1=false;
            }
            if(flag2 && it1.GeSize()>0){
              FillHistogram("pid_De100g",3e3,0,3e4,charge, 4e3,0,8e3,fbeta->fImplant.fPIN1E);
              flag2=false;
            }
          } 
        }
      }
      //}
      //=================================Decay information=========================================//
      bool decayflag = true;
      bool decayflag1 = true;
      int count = 0; // count for decaytime;
      for(auto &it1: fbeta->fDecay){
        count++;
        int i=-1;// count for ggmatrix;
        for(auto &it2: it1.fGe){
          int j1, j2=-1;// count for ggmatrix;
          i++;
          double dt = it1.GetTimestamp() - it2.GetTimestamp();

          for(size_t m=0;m<cut1.size();m++){
            if(cut1[m]->IsInside(charge,fbeta->fImplant.fPIN1E)){
              double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
              dtime = fabs(dtime)/1e6; //unit: ms
              FillHistogram("decaytime_024",cut1[m], 1000,0,500,dtime);
              if(count<2){
                FillHistogram("single_dtime1_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
              }
              if(count<3){
                FillHistogram("single_dtime2_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
              }
              if(count<4){
                FillHistogram("single_dtime3_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
              }
              if(count<5){
                FillHistogram("single_dtime4_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
              }
              if(count<6){
                FillHistogram("single_dtime5_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
              }
              FillHistogram("single_dtime_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            

              if(!it1.GeSize()){
                if(count<2){
                  FillHistogram("decaytime1_024",cut1[m], 1000,0,500,dtime);            
                }
                if(count<3){
                  FillHistogram("decaytime2_024",cut1[m], 1000,0,500,dtime);            
                }
                if(count<4){
                  FillHistogram("decaytime3_024",cut1[m], 1000,0,500,dtime);            
                }
                if(count<5){
                  FillHistogram("decaytime4_024",cut1[m], 1000,0,500,dtime);            
                }
                if(count<6){
                  FillHistogram("decaytime5_024",cut1[m], 1000,0,500,dtime);            
                }
              }

              if(decayflag1 && dtime<=100){
                FillHistogram("decaysize_100",cut1[m], 100,0,100,fbeta->DecaySize());
              }

              //======================================== LHE cut below ==========================================================//
              if(cut2[0]->IsInside(dt,it2.GetEnergy())){
                if(decayflag){FillHistogram("decaysize",cut1[m],100,0,100,fbeta->DecaySize());}
                double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
                dtime = fabs(dtime)/1e6; //unit: ms
                if(count<2){
                  FillHistogram("single_decaytime1_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                }
                if(count<3){
                  //FillHistogram("single_implanttime2_024",cut1[m], 1000,0,500,it1.fImplantTime, 8000,0,4000,it2.GetEnergy());
                  FillHistogram("single_decaytime2_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                  //FillHistogram("decaytime2_024",cut1[m], 1000,0,500,dtime);            
                }
                if(count<4){
                  FillHistogram("single_decaytime3_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                }
                if(count<5){
                  FillHistogram("single_decaytime4_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                }
                if(count<6){
                  FillHistogram("single_decaytime5_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                }
                //FillHistogram("single_implanttime_024",cut1[m], 1000,0,500,it1.fImplantTime, 8000,0,4000,it2.GetEnergy());
                FillHistogram("single_decaytime_024",cut1[m], 1000,0,500,dtime, 8000,0,4000,it2.GetEnergy());            
                //for(auto it3: it1.fGe){
                //  j1++;
                //  if(i<j1){
                //    FillHistogram("ggmat_024",cut1[m], 4000,0,4000,it2.GetEnergy(), 4000,0,4000,it3.GetEnergy());
                //    FillHistogram("ggmat_024",cut1[m], 4000,0,4000,it3.GetEnergy(), 4000,0,4000,it2.GetEnergy());
                //  }
                //}
              }
            }

          } 
          decayflag = false;
          decayflag1 = false;
        }
      }
    }
    if((x%200000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  //SaveHistograms(Form("beta_output%s.root",runnum.c_str())); 
  SaveHistograms(Form("decayNe_output%s.root",runnum.c_str())); 

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
























//================================ Na32 sort ================================//
void Histogram::NaSort(){

  std::string runnum = GetRunNumber(gROOT->GetListOfFiles()->At(0)->GetName());
  runnum = runnum.substr(0,4);

  Beta *fbeta = new Beta;
  gChain->SetBranchAddress("Beta", &fbeta);
  TChannel::ReadDetMapFile();
  long x = 0;
  long n = gChain->GetEntries();

  std::vector<TCutG *> cut1;
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/decaytimecuts.root"); 
  TIter keys1 (mycut1->GetListOfKeys());
  while(TKey *key1 = (TKey*)keys1.Next()) {
    cut1.push_back((TCutG*)key1->ReadObj());
  }

  std::vector<TCutG *> cut2;
  TFile *mycut2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cuts/PID_Na32_cuts.root"); // cuts from pid only Na32 
  TIter keys2 (mycut2->GetListOfKeys());
  while(TKey *key2 = (TKey*)keys2.Next()) {
    cut2.push_back((TCutG*)key2->ReadObj());
  }
  //================== Write new tree in new file ======================//
  TFile *newf = new TFile(Form("Na%s.root",runnum.c_str()),"recreate");
  TTree *tree = new TTree("beta","beta tree");
  Beta *fnewbeta = new Beta;
  tree->Branch("Beta", &fnewbeta); 
  //=====================================================================//


  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    if(fbeta->fImplant.Stopped()){
      for(auto &it1 : fbeta->fDecay){
        int i=-1;//count for ggmat
        for(auto &it2 : it1.fGe){
          int j=-1;//count for ggmat
          i++;
          double dt = it1.GetTimestamp() - it2.GetTimestamp();
          //FillHistogram("decay_gamma_time",1000,-1000,1000,dt, 8000,0,4000,it2.GetEnergy());
          if(cut2[0]->IsInside(dt,it2.GetEnergy())){ // LHE + blob1
            fnewbeta->Set(fbeta->fImplant, fbeta->fDecay);
            tree->Fill();
            fnewbeta->Clear();
            //double dtime = fbeta->fImplant.fDSSDFront[0].GetTimestamp() - it1.fDSSDFront[0].GetTimestamp();
            //FillHistogram("decaytime",1000,0,500,it1.fImplantTime, 8000,0,4000,it2.GetEnergy());
            //FillHistogram("single_implanttime",1000,0,500,it1.fImplantTime, 8000,0,4000,it2.GetEnergy());
            //for(auto &it3 : it1.fGe){
            //  j++; 
            //  if(i<j){
            //    FillHistogram("gg_matrix",8000,0,4000,it2.GetEnergy(), 8000,0,4000,it3.GetEnergy());
            //    FillHistogram("gg_matrix",8000,0,4000,it3.GetEnergy(), 8000,0,4000,it2.GetEnergy());
            //  }
            //}
          }
        }
      }
    }
    if((x%20000)==0){
      printf("on entry %lu / %lu  \r",x,n);
      fflush(stdout);
    }

  }

  printf("   on entry %lu / %lu   \n",x,n);
  SaveHistograms(Form("Na%s.root",runnum.c_str()));
  //tree->Write();
  //newf->Close();


  return;
}



//===========================================================================//





