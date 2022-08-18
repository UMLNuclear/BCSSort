




void beta_tof_check(){


  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Na = (TCutG *)cutf->Get("Na");
  TCutG *Ne = (TCutG *)cutf->Get("Ne");

  TH2D *h1 = new TH2D("tof_Na_before","tof_Na_before",4000,0,4160000,500,11000,16000);
  TH2D *h2 = new TH2D("tof_Na_offsetNa","tof_Na_offsetNa",4000,0,4160000,500,11000,16000);
  TH2D *h4 = new TH2D("tof_Na_offset","tof_Na_offset",4000,0,4160000,500,11000,16000);
  
  TH2D *h5 = new TH2D("tof_Ne_before","tof_Ne_before",4000,0,4160000,500,11000,16000);
  TH2D *h6 = new TH2D("tof_Ne_offsetNa","tof_Ne_offsetNa",4000,0,4160000,500,11000,16000);
  TH2D *h8 = new TH2D("tof_Ne_offset","tof_Ne_offset",4000,0,4160000,500,11000,16000);
  std::vector<double> tof0;  
  std::vector<double> tofNa;  
  std::vector<double> tofNe;  

  long n0 = 0;
  long n = 0;
  long x = 0;
  for(int run=1023;run<1140;run++){
    tof0  = TOFCorrection::Get()->ReadFile(run,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    tofNa = TOFCorrection::Get()->ReadFile(run,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_Na_offset.txt");
    if(tof0.empty()) continue;
    TChain *chan = new TChain("beta");
    chan->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt/beta%i*.root",run));
    Beta *fbeta = new Beta;
    chan->SetBranchAddress("Beta", &fbeta);
    n = chan->GetEntries();
    for(x=0;x<n;x++){
      chan->GetEntry(x);
      double tof = fbeta->fImplant.fI2S;
      double offset   = tof0[0]  + tof;
      double offsetNa = tofNa[0] + tof;
      if(Na->IsInside(fbeta->fImplant.fI2S, fbeta->fImplant.fPIN1E)){
        h1->Fill(x+n0,tof);
        h2->Fill(x+n0,offsetNa);
        h4->Fill(x+n0,offset);
      }
      if(Ne->IsInside(fbeta->fImplant.fI2S, fbeta->fImplant.fPIN1E)){
        h5->Fill(x+n0,tof);
        h6->Fill(x+n0,offsetNa);
        h8->Fill(x+n0,offset);
      }
    }
    n0 += n;
  }
  TFile *newf = new TFile("beta_tof.root","recreate");
  h1->Write();
  h2->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h8->Write();
  newf->Close();

}







