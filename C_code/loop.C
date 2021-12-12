



{ 
  TChannel::ReadDetMapFile(); 
  TChain *chdec = new TChain("decay");
  Decay *fdec = new Decay;
  chdec->Add("run0048/10us/no160175_Esssd600/dec_prompt/decay0048-0*");  
  chdec->SetBranchAddress("Decay", &fdec);

  TChain *chimp = new TChain("implant");
  Implant *fimp = new Implant;
  chimp->Add("run0048/10us/no160175_Esssd600/imp_good/implant0048-0*");  
  chimp->SetBranchAddress("Implant", &fimp);

  TChain *chprob = new TChain("beta");
  TChain *chdelb = new TChain("beta");
  Beta *fprob = new Beta;
  Beta *fdelb = new Beta;
  chprob->Add("run0048/10us/no160175_Esssd600/beta_prompt0048.root");
  chdelb->Add("run0048/10us/no160175_Esssd600/beta_delay0048.root");
  chprob->SetBranchAddress("Beta", &fprob);
  chdelb->SetBranchAddress("Beta", &fdelb);

  long x=0;
  long nimp=chimp->GetEntries();
  long ndec=chdec->GetEntries();
  long nprob=chprob->GetEntries();
  long ndelb=chdelb->GetEntries();

  double low, high=0;
  bool cflag = false;
  int c1 = 0;
  int c2 = 0;
  int c3 = 0;
  int c4 = 0;
  int c5 = 0;
  int c6 = 0;
  int c7 = 0;
  int c8 = 0;
  int c9 = 0;
  int c10 = 0;
  int c11 = 0;
  int c12 = 0;
  int c13 = 0;
  int c14 = 0;
  int c15 = 0;
  int c16 = 0;
  int c17 = 0;
  int c18 = 0;
  int c19 = 0;
  int c20 = 0;
  int c21 = 0;
  int c22 = 0;
  int c23 = 0;
  int c24 = 0;
  int c25 = 0;
  int c26 = 0;
  int c27 = 0;
  TH2D *h1 = new TH2D("single_uncal","single_uncal",16e3,0,32e3, 100,200,300);
  TH2D *h2 = new TH2D("single_cal","single_cal",5e3,0,10e3, 100,200,300);
  TH2D *h3 = new TH2D("single_un_gated","single_un_gated",16e3,0,32e3, 100,200,300);
  //TH1D *h4 = new TH1D("single_cal","single_cal",5e3,0,10e3);

  double dt;
  for(x=0;x<nprob;x++){
    chprob->GetEntry(x);
    bool flag = false;
    if(fprob->DecaySize()>0) c1++;
    for(auto &it1:fprob->fDecay){
      if(it1.fGe.size()==0) continue;
      for(auto &it2:it1.fGe){
        if(it2.GetEnergy()<10 || it2.GetEnergy()<4000) continue;
        flag = true;
        break;
      }
    }
    if(flag) c2++;
  }
  for(x=0;x<ndelb;x++){
    chdelb->GetEntry(x);
    bool flag = false;
    if(fdelb->DecaySize()==0) continue;
    c3++;
    for(auto &it1:fdelb->fDecay){
      if(it1.fGe.size()==0) continue;
      for(auto &it2:it1.fGe){
        if(it2.GetEnergy()<10 || it2.GetEnergy()<4000) continue;
        flag = true;
        break;
      }
    }
    if(flag) c4++;
  }
  
  std::cout<< "total entries of prompt correlation = "<< nprob<< std::endl;
  std::cout<< "Prompt correlation with nonzero decay = " << c1 << std::endl;
  std::cout<< "beta count with decays which have good gamma = " << c2 << std::endl;
  std::cout<<std::endl;
  std::cout<< "total entries of delay correlation = "<< ndelb<< std::endl;
  std::cout<< "Delay correlation with nonzero decay = " << c3 << std::endl;
  std::cout<< "beta count with decays which have good gamma = " << c4 << std::endl;
} 
