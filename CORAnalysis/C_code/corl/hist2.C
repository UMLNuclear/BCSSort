



void hist2(){

  TChannel::ReadDetMapFile();
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");

  TH2D *h0 = 0;  
  TH2D *h1[6];  
  TH2D *e0 = 0;  
  TH2D *e1[6];  
  TH2D *te1[6];
  TH1D *t1[6];
  TH1D *g0 = 0;
  TH2D *p0[6];
  TH2D *pb[6];
  
  double game[6] = {546, 885, 2244, 150, 1516, 2063};
  double gam1[6] = {544, 883, 2242, 147, 1514, 2061};
  double gam2[6] = {549, 888, 2247, 152, 1519, 2066};
  double gbg1[6] = {549, 888, 2247, 142, 1519, 2066};
  double gbg2[6] = {554, 893, 2252, 147, 1524, 2071};

  for(int i=0;i<6;i++){
    h1[i]  = new TH2D(Form("pid%i",i),Form("pid %s",blob[i]->GetName()), 1000,8000,18000,300,4200,7200); 
    e1[i]  = new TH2D(Form("e%i",i)  ,Form("dE(LGF) vs E(SSSD) %s",blob[i]->GetName()), 400,0,32000,400,0,4000); 
    te1[i] = new TH2D(Form("te%i",i)  ,Form("gammaE(keV) vs decaytime(ms) %s",blob[i]->GetName()), 100,0,100,4000,0,4000); 
    t1[i]  = new TH1D(Form("t%i",i)  ,Form("decaytime(ms) %s",blob[i]->GetName()), 10000,0,1000); 
    p0[i]  = new TH2D(Form("pidg%i",i),Form("pid gamma gated %.2fkeV",game[i]), 1000,8000,18000,300,4200,7200); 
    pb[i]  = new TH2D(Form("pidb%i",i),Form("pid(bg) gamma gated %.2fkeV",game[i]), 1000,8000,18000,300,4200,7200); 
  }    

  for(int runnum = 1040;runnum<1050;runnum++){
    std::vector<double> tofpar= TOFCorrection::Get()->ReadFile(runnum);
    std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);

    if(!h0){ h0 = new TH2D("pid"     , "pid"     , 1000,8000,18000,300,4200,7200); }
    if(!e0){ e0 = new TH2D("e",      "dE(LGF) vs E(SSSD)"     , 400,0,32000,400,0,4000); }
    if(!g0){ g0 = new TH1D("decaysize", "decaysize", 1000,0,1000);}


    long nentries = chan->GetEntries();
    long x = 0;
    double pin1t, pin1e, tof;
    int decaysize = 0;
    bool flag_blob[6];

    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0){
        if(fevent->LGFSize()>10) continue;
        for(int i=0;i<6;i++) { flag_blob[i] = false; }
        if(x!=0){ g0->Fill(decaysize); }
        decaysize = 0;
        pin1e = fevent->Pin1E();
        pin1t = fevent->Pin1T();
        tof   = fevent->I2S();
        tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
        h0->Fill(tof, pin1e);
        for(int i=0;i<6;i++){
          if(blob[i]->IsInside(tof,pin1e)){
            h1[i]->Fill(tof,pin1e);
            flag_blob[i] = true;
          }
        }
        if(fevent->SSSDSize()>0){
          for(int j=0;j<fevent->SSSDSize();j++){
            int    sssdn = fevent->SSSD()[j].GetStrip()-1;
            double sssde = fevent->SSSD()[j].GetCharge();
            if(!cutstrip->IsInside(sssdn, (fevent->LGFMax().GetStrip()-1))) continue;
            e0->Fill(sssde, fevent->LGFMax().GetEnergy());
            for(int i=0;i<6;i++){
              if(flag_blob[i]) { e1[i]->Fill(sssde, fevent->LGFMax().GetEnergy()); } 
            }  
          }
        }
      }else{
        decaysize++;
        double dt = fevent->HGF()[0].GetTimestamp() - pin1t;
        dt = dt/1000000.;
        for(int i=0;i<6;i++){
          if(flag_blob[i]) { t1[i]->Fill(dt); }
        }
        if(fevent->HPGeSize()>0){
          for(int i=0;i<fevent->HPGeSize();i++){
            double e = fevent->HPGe()[i].GetEnergy();
            if(e<10 || e>4000) continue;
            for(int i=0;i<6;i++){
              if(flag_blob[i]) { te1[i]->Fill(dt, e); }
              if(dt>30) continue;
              if(e>=gam1[i] && e<=gam2[i]){ p0[i]->Fill(tof, pin1e); }
              if(e>=gbg1[i] && e<=gbg2[i]){ pb[i]->Fill(tof, pin1e); }
            }
          }
        } 
      }


      if((x%5000)==0){
        printf("on entry %lu / %lu \r",x,nentries);
        fflush(stdout);
      }
    }


    printf("on entry %lu / %lu \n",x,nentries);
  }

  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/cor2_hist.root","recreate");
  h0->Write();
  e0->Write();
  g0->Write();
  for(int i=0;i<6;i++){
    h1[i] ->Write();
    e1[i] ->Write();
    te1[i]->Write();
    t1[i] ->Write();
    p0[i] ->Write();
    pb[i] ->Write();
  }
  newf->Close();

  

}
