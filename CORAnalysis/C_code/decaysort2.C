

void hist(int runnum=1040){
  TChannel::ReadDetMapFile();
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/hgfc_hbc_cuts.root");
  std::map<std::pair<int,int>, TCutG *> hgc0;
  std::map<std::pair<int,int>, TCutG *> hgc1;
  std::map<std::pair<int,int>, TCutG *> hgc2;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      std::pair<int,int> temp = std::make_pair(i+1,j+1);
      hgc0[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,0));
      hgc1[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,1));
      hgc2[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,2));
    }
  }

  TH2D *cc0[9]; // no time gate;
  TH2D *cc1[9]; // with time gate;
  TH2D *ss[3];  // hgf.size vs hgb.size with gates;
  TH2D *hm[3];  // hitmap dt gated
  TH2D *gt[3];  // gamma vs dt = gamma-hgfmaxt, dt gated

  for(int i=0;i<3;i++){
    ss[i] = new TH2D(Form("ss%i",i),Form("hgf.size vs hgb.size() only with pixels in time gate for c.zone%i",i),40,0,40,40,0,40);
    hm[i] = new TH2D(Form("hm%i",i),Form("hitmap time gated in c.zone%i",i), 40,0,40,40,0,40);
    gt[i] = new TH2D(Form("gt%i",i),Form("dt(ns)=gamma-hgfmaxt vs gammaE(keV) with dt=[30,140)ns in c.zone%i",i),1000,-5000,5000,4000,0,4000);
    for(int j=0;j<3;j++){
      cc0[3*i+j] = new TH2D(Form("cc0_%i%i",i+1,j+1),Form("hgf(%i).max vs hgb(%i).max",i+1,j+1),500,0,25000,500,0,25000);
      cc1[3*i+j] = new TH2D(Form("cc1_%i%i",i+1,j+1),Form("hgf(%i).max vs hgb(%i).max with time gate",i+1,j+1),500,0,25000,500,0,25000);
    }
  }

  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/decay/highgain/event%i*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long x=0;
  long nentries = chan->GetEntries();

  std::vector<int> vf;
  std::vector<int> vb;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->HGFMax().GetEnergy()>0 && fevent->HGBMax().GetEnergy()>0){
      int hgfmaxn = fevent->HGFMax().GetStrip()-1;
      int hgbmaxn = fevent->HGBMax().GetStrip()-1;
      double hgfmaxc = fevent->HGFMax().GetCharge();
      double hgbmaxc = fevent->HGBMax().GetCharge();
      double hgfmaxt = fevent->HGFMax().GetTimestamp();
      double hgbmaxt = fevent->HGBMax().GetTimestamp();
      int car = 3*(hgfmaxn/16)+(hgbmaxn/16);
      cc0[car]->Fill(hgfmaxc, hgbmaxc);
      double dtmax = hgfmaxt - hgbmaxt;
      if(dtmax>=30 && dtmax<140){
        cc1[car]->Fill(hgfmaxc, hgbmaxc);
        std::pair<int,int> temp = std::make_pair((hgfmaxn/16)+1, (hgbmaxn/16)+1);
        if(hgc0[temp]->IsInside(hgfmaxc, hgbmaxc)){ ss[0]->Fill(fevent->HGFSize(10,25000), fevent->HGBSize(10,25000)); hm[0]->Fill(hgfmaxn,hgbmaxn);}
        if(hgc1[temp]->IsInside(hgfmaxc, hgbmaxc)){ ss[1]->Fill(fevent->HGFSize(10,25000), fevent->HGBSize(10,25000)); hm[1]->Fill(hgfmaxn,hgbmaxn);}
        if(hgc2[temp]->IsInside(hgfmaxc, hgbmaxc)){ ss[2]->Fill(fevent->HGFSize(10,25000), fevent->HGBSize(10,25000)); hm[2]->Fill(hgfmaxn,hgbmaxn);}
        for(int m=0;m<fevent->HPGeSize();m++){
          double e = fevent->HPGe()[m].GetEnergy();
          double decaytime = fevent->HPGe()[m].GetTimestamp() - hgfmaxt;
          if(e<10 || e>4000) continue;
          if(hgc0[temp]->IsInside(hgfmaxc, hgbmaxc)){ gt[0]->Fill(decaytime,e);}
          if(hgc1[temp]->IsInside(hgfmaxc, hgbmaxc)){ gt[1]->Fill(decaytime,e);}
          if(hgc2[temp]->IsInside(hgfmaxc, hgbmaxc)){ gt[2]->Fill(decaytime,e);}
        }// gamma loop end
      }// time gated
    }// hgfmax && hgbmax exists
    


    if((x%5000)==0){
      printf("run%i.root: on entry %lu / %lu \r", runnum,x,nentries);
      fflush(stdout);
    }
  }
  printf("run%i.root: on entry %lu / %lu \n", runnum,x,nentries);
  TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/highgain/hgfc_hgbc1/dechist%i.root",runnum),"recreate");
  for(int i=0;i<3;i++){
    ss[i] -> Write();
    hm[i] -> Write();
    gt[i] -> Write();
    for(int j=0;j<3;j++){
      cc0[3*i+j] -> Write();
      cc1[3*i+j] -> Write();
    }
  }
  newf->Close();

}

void hist_gamgated(int runnum=1040){

  TChannel::ReadDetMapFile();
  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/decay/highgain/event%i*.root",runnum));

  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long x = 0;
  long nentries = chan->GetEntries();
  
  TH2D *h0[9]; // no gate;
  TH2D *h10[9]; // dt gate;
  TH2D *h1[9]; // have gamma;
  TH2D *h2[9]; // gamma gate at 885; (882,886)keV
  TH2D *h3[9]; // gamma gate at 885; (882,886)keV inside dt gate;
  TH2D *h4[9]; // gamma gate at 885BG; (900,906)keV;
  TH2D *h5[9]; // gamma gate at 885BG; (900,906)keV inside dt gate;
  TH2D *h6[9]; // gamma gate at 150; (146,152)keV
  TH2D *h7[9]; // gamma gate at 150; (146,152)keV inside dt gate;
  TH2D *h8[9]; // gamma gate at 150BG; (155,161)keV;
  TH2D *h9[9]; // gamma gate at 150BG; (155,161)keV inside dt gate;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      h0[3*i+j] = new TH2D(Form("c%i%i",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i)",i,j), 500,0,25000,500,0,25000);
      h10[3*i+j] = new TH2D(Form("ct%i%i",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) within dt=[30,140)ns",i,j) , 500,0,25000,500,0,25000);
      h1[3*i+j] = new TH2D(Form("cg%i%i",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) has gamma",i,j)    , 500,0,25000,500,0,25000);
      h2[3*i+j] = new TH2D(Form("cg%i%i_885",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 885",i,j), 500,0,25000,500,0,25000);
      h3[3*i+j] = new TH2D(Form("cg%i%i_885dt",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 885 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      h4[3*i+j] = new TH2D(Form("cg%i%i_885bg",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG",i,j), 500,0,25000,500,0,25000);
      h5[3*i+j] = new TH2D(Form("cg%i%i_885dtbg",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      h6[3*i+j] = new TH2D(Form("cg%i%i_150",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 150",i,j), 500,0,25000,500,0,25000);
      h7[3*i+j] = new TH2D(Form("cg%i%i_150dt",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 150 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      h8[3*i+j] = new TH2D(Form("cg%i%i_150bg",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG",i,j), 500,0,25000,500,0,25000);
      h9[3*i+j] = new TH2D(Form("cg%i%i_150dtbg",i+1,j+1), Form("hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
    }
  }
  
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->HGFSize(10,25000)==0 || fevent->HGBSize(10,25000)==0) continue;
    int    hgfn = fevent->HGFMax().GetStrip()-1;
    double hgft = fevent->HGFMax().GetTimestamp();
    double hgfc = fevent->HGFMax().GetCharge();
    if(hgfc<10 || hgfc>25000) continue;
    int    hgbn = fevent->HGBMax().GetStrip()-1;
    double hgbt = fevent->HGBMax().GetTimestamp();
    double hgbc = fevent->HGBMax().GetCharge();
    if(hgbc<10 || hgbc>25000) continue;
    int car = 3*(hgfn/16)+(hgbn/16);
    h0[car]->Fill(hgfc, hgbc);
    double dt = hgft - hgbt;
    if(dt>=30 && dt<140){
      h10[car]->Fill(hgfc,hgbc);
    }
    if(fevent->HPGeSize()>0){
      h1[car]->Fill(hgfc, hgbc);
      for(int m=0;m<fevent->HPGeSize();m++){
        double e = fevent->HPGe()[m].GetEnergy();
        if(e>882 && e<888){ // gamma gated 885keV
          h2[car]->Fill(hgfc, hgbc);
          if(dt>=30 && dt<140){
            h3[car]->Fill(hgfc, hgbc);
          }
        }
        if(e>900 && e<906){ // gamma gated 885keVBG
          h4[car]->Fill(hgfc, hgbc);
          if(dt>=30 && dt<140){
            h5[car]->Fill(hgfc, hgbc);
          }
        }
        if(e>146 && e<152){ // gamma gated 150keV
          h6[car]->Fill(hgfc, hgbc);
          if(dt>=30 && dt<140){
            h7[car]->Fill(hgfc, hgbc);
          }
        }
        if(e>155 && e<161){ // gamma gated 150keVBG
          h8[car]->Fill(hgfc, hgbc);
          if(dt>=30 && dt<140){
            h9[car]->Fill(hgfc, hgbc);
          }
        }
      }//gamma loop end
    }// if gamma exist
    if((x%20000)==0){
      printf("run%i.root: %lu / %lu \r", runnum, x,nentries);
      fflush(stdout);
    }
  }
    
  printf("run%i.root: %lu / %lu \n", runnum, x,nentries);
  TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/highgain/hgfc_hgbc/maxc/dechist%i.root", runnum), "recreate"); 
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      h0[3*i+j] -> Write();
      h1[3*i+j] -> Write();
      h2[3*i+j] -> Write();
      h3[3*i+j] -> Write();
      h4[3*i+j] -> Write();
      h5[3*i+j] -> Write();
      h6[3*i+j] -> Write();
      h7[3*i+j] -> Write();
      h8[3*i+j] -> Write();
      h9[3*i+j] -> Write();
      h10[3*i+j] -> Write();
    }
  }
  newf->Close();

}


void hist_gamgated_IsotopeGated(int runnum=1040){

  TChannel::ReadDetMapFile();

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/hgfc_hbc_cuts.root");
  std::map<std::pair<int,int>, TCutG *> hgc0;
  std::map<std::pair<int,int>, TCutG *> hgc1;
  std::map<std::pair<int,int>, TCutG *> hgc2;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      std::pair<int,int> temp = std::make_pair(i+1,j+1);
      hgc0[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,0));
      hgc1[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,1));
      hgc2[temp] = (TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,2));
    }
  }
  
  TChain *chan = new TChain("beta");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/correlation4/beta%i*.root",runnum));
  Beta *fbeta = new Beta;
  chan->SetBranchAddress("Beta", &fbeta);
  long x = 0;
  long nentries = chan->GetEntries();
  
  TH2D *a0[9]; // no gate;
  TH2D *a10[9]; // dt gate;
  TH2D *a1[9]; // have gamma;
  TH2D *a2[9]; // gamma gate at 885; (882,886)keV
  TH2D *a3[9]; // gamma gate at 885; (882,886)keV inside dt gate;
  TH2D *a4[9]; // gamma gate at 885BG; (900,906)keV;
  TH2D *a5[9]; // gamma gate at 885BG; (900,906)keV inside dt gate;
  TH2D *a6[9]; // gamma gate at 150; (146,152)keV
  TH2D *a7[9]; // gamma gate at 150; (146,152)keV inside dt gate;
  TH2D *a8[9]; // gamma gate at 150BG; (155,161)keV;
  TH2D *a9[9]; // gamma gate at 150BG; (155,161)keV inside dt gate;
  TH2D *a11[3]; // hitmap in different hg.c zone;
  TH2D *a12[3]; // dt = gamma.t - hgfmax.t vs gamma.E in different hg.c zones;
  TH2D *e0[9]; // no gate;
  TH2D *e10[9]; // dt gate;
  TH2D *e1[9]; // have gamma;
  TH2D *e2[9]; // gamma gate at 885; (882,886)keV
  TH2D *e3[9]; // gamma gate at 885; (882,886)keV inside dt gate;
  TH2D *e4[9]; // gamma gate at 885BG; (900,906)keV;
  TH2D *e5[9]; // gamma gate at 885BG; (900,906)keV inside dt gate;
  TH2D *e6[9]; // gamma gate at 150; (146,152)keV
  TH2D *e7[9]; // gamma gate at 150; (146,152)keV inside dt gate;
  TH2D *e8[9]; // gamma gate at 150BG; (155,161)keV;
  TH2D *e9[9]; // gamma gate at 150BG; (155,161)keV inside dt gate;
  TH2D *e11[3]; // hitmap in different hg.c zone;
  TH2D *e12[3]; // dt = gamma.t - hgfmax.t vs gamma.E in different hg.c zones;
  for(int i=0;i<3;i++){
    a11[i] = new TH2D(Form("ahm%i",i),Form("32Na: hitmap time gated in c.zone%i",i), 40,0,40,40,0,40);
    a12[i] = new TH2D(Form("agt%i",i),Form("32Na: dt(ns)=gamma-hgfmaxt vs gammaE(keV) with dt=[30,140)ns in c.zone%i",i),1000,-5000,5000,4000,0,4000);
    e11[i] = new TH2D(Form("ehm%i",i),Form("Ne: hitmap time gated in c.zone%i",i), 40,0,40,40,0,40);
    e12[i] = new TH2D(Form("egt%i",i),Form("Ne: dt(ns)=gamma-hgfmaxt vs gammaE(keV) with dt=[30,140)ns in c.zone%i",i),1000,-5000,5000,4000,0,4000);
    for(int j=0;j<3;j++){
      a0[3*i+j] = new TH2D(Form("a%i%i",i+1,j+1)         , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i)",i,j), 500,0,25000,500,0,25000);
      a10[3*i+j]= new TH2D(Form("at%i%i",i+1,j+1)        , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) within dt=[30,140)ns",i,j) , 500,0,25000,500,0,25000);
      a1[3*i+j] = new TH2D(Form("ag%i%i",i+1,j+1)        , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) has gamma",i,j)    , 500,0,25000,500,0,25000);
      a2[3*i+j] = new TH2D(Form("ag%i%i_885",i+1,j+1)    , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 885",i,j), 500,0,25000,500,0,25000);
      a3[3*i+j] = new TH2D(Form("ag%i%i_885dt",i+1,j+1)  , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 885 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      a4[3*i+j] = new TH2D(Form("ag%i%i_885bg",i+1,j+1)  , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG",i,j), 500,0,25000,500,0,25000);
      a5[3*i+j] = new TH2D(Form("ag%i%i_885dtbg",i+1,j+1), Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      a6[3*i+j] = new TH2D(Form("ag%i%i_150",i+1,j+1)    , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 150",i,j), 500,0,25000,500,0,25000);
      a7[3*i+j] = new TH2D(Form("ag%i%i_150dt",i+1,j+1)  , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 150 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      a8[3*i+j] = new TH2D(Form("ag%i%i_150bg",i+1,j+1)  , Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG",i,j), 500,0,25000,500,0,25000);
      a9[3*i+j] = new TH2D(Form("ag%i%i_150dtbg",i+1,j+1), Form("32Na: hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      e0[3*i+j] = new TH2D(Form("e%i%i",i+1,j+1)         , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i)",i,j), 500,0,25000,500,0,25000);
      e10[3*i+j]= new TH2D(Form("et%i%i",i+1,j+1)        , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) within dt=[30,140)ns",i,j) , 500,0,25000,500,0,25000);
      e1[3*i+j] = new TH2D(Form("eg%i%i",i+1,j+1)        , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) has gamma",i,j)    , 500,0,25000,500,0,25000);
      e2[3*i+j] = new TH2D(Form("eg%i%i_885",i+1,j+1)    , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 885",i,j), 500,0,25000,500,0,25000);
      e3[3*i+j] = new TH2D(Form("eg%i%i_885dt",i+1,j+1)  , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 885 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      e4[3*i+j] = new TH2D(Form("eg%i%i_885bg",i+1,j+1)  , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG",i,j), 500,0,25000,500,0,25000);
      e5[3*i+j] = new TH2D(Form("eg%i%i_885dtbg",i+1,j+1), Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 885BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      e6[3*i+j] = new TH2D(Form("eg%i%i_150",i+1,j+1)    , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 150",i,j), 500,0,25000,500,0,25000);
      e7[3*i+j] = new TH2D(Form("eg%i%i_150dt",i+1,j+1)  , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 150 within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
      e8[3*i+j] = new TH2D(Form("eg%i%i_150bg",i+1,j+1)  , Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG",i,j), 500,0,25000,500,0,25000);
      e9[3*i+j] = new TH2D(Form("eg%i%i_150dtbg",i+1,j+1), Form("Ne: hgfmax.c(%i) vs hgbmax.c(%i) gated 150BG within dt=[30,140)ns",i,j), 500,0,25000,500,0,25000);
    }
  }
  
  for(x=0;x<nentries;x++){
    fbeta->Clear();
    chan->GetEntry(x);
    double tof   = fbeta->fImplant.fI2S;
    double pin1e = fbeta->fImplant.fPIN1E;
    if(blob[1]->IsInside(tof, pin1e)){
      if(fbeta->DecaySize()==0) continue;
      for(int d=0;d<fbeta->DecaySize();d++){
        Decay fevent;
        fevent.Clear();
        fevent = fbeta->fDecay[d]; 
        if(fevent.FrontSize(10,25000)==0 || fevent.FrontSize(10,25000)==0) continue;
        int    hgfn = fevent.HGFMax().GetStrip()-1;
        double hgft = fevent.HGFMax().GetTimestamp();
        double hgfc = fevent.HGFMax().GetCharge();
        if(hgfc<10 || hgfc>25000) continue;
        int    hgbn = fevent.HGBMax().GetStrip()-1;
        double hgbt = fevent.HGBMax().GetTimestamp();
        double hgbc = fevent.HGBMax().GetCharge();
        if(hgbc<10 || hgbc>25000) continue;
        int car = 3*(hgfn/16)+(hgbn/16);
        a0[car]->Fill(hgfc, hgbc);
        std::pair<int,int> temp = std::make_pair((hgfn/16)+1, (hgbn/16)+1);
        if(hgc0[temp]->IsInside(hgfc, hgbc)){ a11[0]->Fill(hgfn,hgbn);}
        if(hgc1[temp]->IsInside(hgfc, hgbc)){ a11[1]->Fill(hgfn,hgbn);}
        if(hgc2[temp]->IsInside(hgfc, hgbc)){ a11[2]->Fill(hgfn,hgbn);}
        double dt = hgft - hgbt;
        if(dt>=30 && dt<140){
          a10[car]->Fill(hgfc,hgbc);
        }
        if(fevent.GeSize()>0){
          a1[car]->Fill(hgfc, hgbc);
          for(int m=0;m<fevent.GeSize();m++){
            double e = fevent.fGe[m].GetEnergy();
            double decayt = fevent.fGe[m].GetTimestamp() - hgft;
            if(e<10 || e>4000) continue;
            if(hgc0[temp]->IsInside(hgfc, hgbc)){ a12[0]->Fill(decayt,e);}
            if(hgc1[temp]->IsInside(hgfc, hgbc)){ a12[1]->Fill(decayt,e);}
            if(hgc2[temp]->IsInside(hgfc, hgbc)){ a12[2]->Fill(decayt,e);}
            if(e>882 && e<888){ // gamma gated 885keV
              a2[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                a3[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>900 && e<906){ // gamma gated 885keVBG
              a4[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                a5[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>146 && e<152){ // gamma gated 150keV
              a6[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                a7[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>155 && e<161){ // gamma gated 150keVBG
              a8[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                a9[car]->Fill(hgfc, hgbc);
              }
            }
          }//gamma loop end
        }// if gamma exist
      }// decay loop end
    }// 32Na implant gate end
    if(blob[3]->IsInside(tof, pin1e)){
      if(fbeta->DecaySize()==0) continue;
      for(int d=0;d<fbeta->DecaySize();d++){
        Decay fevent;
        fevent.Clear();
        fevent = fbeta->fDecay[d]; 
        if(fevent.FrontSize(10,25000)==0 || fevent.FrontSize(10,25000)==0) continue;
        int    hgfn = fevent.HGFMax().GetStrip()-1;
        double hgft = fevent.HGFMax().GetTimestamp();
        double hgfc = fevent.HGFMax().GetCharge();
        if(hgfc<10 || hgfc>25000) continue;
        int    hgbn = fevent.HGBMax().GetStrip()-1;
        double hgbt = fevent.HGBMax().GetTimestamp();
        double hgbc = fevent.HGBMax().GetCharge();
        if(hgbc<10 || hgbc>25000) continue;
        int car = 3*(hgfn/16)+(hgbn/16);
        e0[car]->Fill(hgfc, hgbc);
        std::pair<int,int> temp = std::make_pair((hgfn/16)+1, (hgbn/16)+1);
        if(hgc0[temp]->IsInside(hgfc, hgbc)){ e11[0]->Fill(hgfn,hgbn);}
        if(hgc1[temp]->IsInside(hgfc, hgbc)){ e11[1]->Fill(hgfn,hgbn);}
        if(hgc2[temp]->IsInside(hgfc, hgbc)){ e11[2]->Fill(hgfn,hgbn);}
        double dt = hgft - hgbt;
        if(dt>=30 && dt<140){
          e10[car]->Fill(hgfc,hgbc);
        }
        if(fevent.GeSize()>0){
          e1[car]->Fill(hgfc, hgbc);
          for(int m=0;m<fevent.GeSize();m++){
            double e = fevent.fGe[m].GetEnergy();
            double decayt = fevent.fGe[m].GetTimestamp() - hgft;
            if(e<10 || e>4000) continue;
            if(hgc0[temp]->IsInside(hgfc, hgbc)){ e12[0]->Fill(decayt,e);}
            if(hgc1[temp]->IsInside(hgfc, hgbc)){ e12[1]->Fill(decayt,e);}
            if(hgc2[temp]->IsInside(hgfc, hgbc)){ e12[2]->Fill(decayt,e);}
            if(e>882 && e<888){ // gamma gated 885keV
              e2[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                e3[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>900 && e<906){ // gamma gated 885keVBG
              e4[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                e5[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>146 && e<152){ // gamma gated 150keV
              e6[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                e7[car]->Fill(hgfc, hgbc);
              }
            }
            if(e>155 && e<161){ // gamma gated 150keVBG
              e8[car]->Fill(hgfc, hgbc);
              if(dt>=30 && dt<140){
                e9[car]->Fill(hgfc, hgbc);
              }
            }
          }//gamma loop end
        }// if gamma exist
      }// decay loop end
    }// 31+30Ne implant gate end
    if((x%20000)==0){
      printf("run%i.root: %lu / %lu \r", runnum, x,nentries);
      fflush(stdout);
    }
  }
    
  printf("run%i.root: %lu / %lu \n", runnum, x,nentries);
  TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/highgain/hgfc_hgbc/maxc/IsotopeGated/cor4hist%i.root", runnum), "recreate"); 
  for(int i=0;i<3;i++){
    a11[i] -> Write();
    a12[i] -> Write();
    e11[i] -> Write();
    e12[i] -> Write();
    for(int j=0;j<3;j++){
      a0[3*i+j] -> Write();
      a1[3*i+j] -> Write();
      a2[3*i+j] -> Write();
      a3[3*i+j] -> Write();
      a4[3*i+j] -> Write();
      a5[3*i+j] -> Write();
      a6[3*i+j] -> Write();
      a7[3*i+j] -> Write();
      a8[3*i+j] -> Write();
      a9[3*i+j] -> Write();
      a10[3*i+j] -> Write();
      e0[3*i+j] -> Write();
      e1[3*i+j] -> Write();
      e2[3*i+j] -> Write();
      e3[3*i+j] -> Write();
      e4[3*i+j] -> Write();
      e5[3*i+j] -> Write();
      e6[3*i+j] -> Write();
      e7[3*i+j] -> Write();
      e8[3*i+j] -> Write();
      e9[3*i+j] -> Write();
      e10[3*i+j] -> Write();
    }
  }
  newf->Close();

}
