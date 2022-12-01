





void decaysort(){


  TChannel::ReadDetMapFile();
  TChain *chan = new TChain("event");
  chan->Add("/home/zhu/packages/BCSSort/CORAnalysis/example/data/decay/highgain/event104*.root");
  //chan->Add("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event104*.root");
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);  

  long nentries = chan->GetEntries();
  long x = 0;

  TH2D *e00 = new TH2D("e00", "hgfe    vs hgbe   ", 400,0,16000,400,0,16000);
  TH2D *e11 = new TH2D("e11", "hgfe(1) vs hgbe(1)", 400,0,16000,400,0,16000);
  TH2D *e12 = new TH2D("e12", "hgfe(1) vs hgbe(2)", 400,0,16000,400,0,16000);
  TH2D *e13 = new TH2D("e13", "hgfe(1) vs hgbe(3)", 400,0,16000,400,0,16000);
  TH2D *e21 = new TH2D("e21", "hgfe(2) vs hgbe(1)", 400,0,16000,400,0,16000);
  TH2D *e22 = new TH2D("e22", "hgfe(2) vs hgbe(2)", 400,0,16000,400,0,16000);
  TH2D *e23 = new TH2D("e23", "hgfe(2) vs hgbe(3)", 400,0,16000,400,0,16000);
  TH2D *e31 = new TH2D("e31", "hgfe(3) vs hgbe(1)", 400,0,16000,400,0,16000);
  TH2D *e32 = new TH2D("e32", "hgfe(3) vs hgbe(2)", 400,0,16000,400,0,16000);
  TH2D *e33 = new TH2D("e33", "hgfe(3) vs hgbe(3)", 400,0,16000,400,0,16000);

  TH1D *t00 = new TH1D("t00", "dt(ns) = hgft    - hgbt   ", 1000,-5000,5000);
  TH1D *t11 = new TH1D("t11", "dt(ns) = hgft(1) - hgbt(1)", 1000,-5000,5000);
  TH1D *t12 = new TH1D("t12", "dt(ns) = hgft(1) - hgbt(2)", 1000,-5000,5000);
  TH1D *t13 = new TH1D("t13", "dt(ns) = hgft(1) - hgbt(3)", 1000,-5000,5000);
  TH1D *t21 = new TH1D("t21", "dt(ns) = hgft(2) - hgbt(1)", 1000,-5000,5000);
  TH1D *t22 = new TH1D("t22", "dt(ns) = hgft(2) - hgbt(2)", 1000,-5000,5000);
  TH1D *t23 = new TH1D("t23", "dt(ns) = hgft(2) - hgbt(3)", 1000,-5000,5000);
  TH1D *t31 = new TH1D("t31", "dt(ns) = hgft(3) - hgbt(1)", 1000,-5000,5000);
  TH1D *t32 = new TH1D("t32", "dt(ns) = hgft(3) - hgbt(2)", 1000,-5000,5000);
  TH1D *t33 = new TH1D("t33", "dt(ns) = hgft(3) - hgbt(3)", 1000,-5000,5000);


  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->Pin1E()<0){
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        for(int i=0;i<fevent->HGFSize();i++){
          bool f1 = false;
          bool f2 = false;
          bool f3 = false;
          int    hgfn = fevent->HGF()[i].GetStrip()-1;
          double hgfe = fevent->HGF()[i].GetCharge();
          double hgft = fevent->HGF()[i].GetTimestamp();
          if(fevent->Range(hgfn,0,15))  f1 = true;
          if(fevent->Range(hgfn,16,31)) f2 = true;
          if(fevent->Range(hgfn,32,39)) f3 = true;
          for(int j=0;j<fevent->HGBSize();j++){
            bool b1 = false;
            bool b2 = false;
            bool b3 = false;
            int    hgbn = fevent->HGB()[j].GetStrip()-1;
            double hgbe = fevent->HGB()[j].GetCharge();
            double hgbt = fevent->HGB()[j].GetTimestamp();
            if(fevent->Range(hgbn,0,15))  b1 = true;
            if(fevent->Range(hgbn,16,31)) b2 = true;
            if(fevent->Range(hgbn,32,39)) b3 = true;
            e00->Fill(hgfe,hgbe);
            t00->Fill((hgft-hgbt)); 
            if(f1 && b1){ e11->Fill(hgfe, hgbe); t11->Fill((hgft-hgbt)); } 
            if(f1 && b2){ e12->Fill(hgfe, hgbe); t12->Fill((hgft-hgbt)); } 
            if(f1 && b3){ e13->Fill(hgfe, hgbe); t13->Fill((hgft-hgbt)); } 
            if(f2 && b1){ e21->Fill(hgfe, hgbe); t21->Fill((hgft-hgbt)); } 
            if(f2 && b2){ e22->Fill(hgfe, hgbe); t22->Fill((hgft-hgbt)); } 
            if(f2 && b3){ e23->Fill(hgfe, hgbe); t23->Fill((hgft-hgbt)); } 
            if(f3 && b1){ e31->Fill(hgfe, hgbe); t31->Fill((hgft-hgbt)); } 
            if(f3 && b2){ e32->Fill(hgfe, hgbe); t32->Fill((hgft-hgbt)); } 
            if(f3 && b3){ e33->Fill(hgfe, hgbe); t33->Fill((hgft-hgbt)); } 

          }
        }
      }

    }

    if((x%5000)==0){
      printf("on entry = %lu / %lu \r",x,nentries);
      fflush(stdout);
    }

  }

  printf("on entry %lu / %lu \n", x, nentries);

  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/no_correlation.root","recreate");
  //TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2.root","recreate");
  e00->Write();
  e11->Write();
  e12->Write();
  e13->Write();
  e21->Write();
  e22->Write();
  e23->Write();
  e31->Write();
  e32->Write();
  e33->Write();
  t00->Write();
  t11->Write();
  t12->Write();
  t13->Write();
  t21->Write();
  t22->Write();
  t23->Write();
  t31->Write();
  t32->Write();
  t33->Write();
  newf->Close();

}




void decaysort_cor(){
  
  TChannel::ReadDetMapFile();

  TH2D *c00 = new TH2D("c00", "hgfc    vs hgbc   ", 400,0,16000,400,0,16000);
  TH2D *c11 = new TH2D("c11", "hgfc(1) vs hgbc(1)", 400,0,16000,400,0,16000);
  TH2D *c12 = new TH2D("c12", "hgfc(1) vs hgbc(2)", 400,0,16000,400,0,16000);
  TH2D *c13 = new TH2D("c13", "hgfc(1) vs hgbc(3)", 400,0,16000,400,0,16000);
  TH2D *c21 = new TH2D("c21", "hgfc(2) vs hgbc(1)", 400,0,16000,400,0,16000);
  TH2D *c22 = new TH2D("c22", "hgfc(2) vs hgbc(2)", 400,0,16000,400,0,16000);
  TH2D *c23 = new TH2D("c23", "hgfc(2) vs hgbc(3)", 400,0,16000,400,0,16000);
  TH2D *c31 = new TH2D("c31", "hgfc(3) vs hgbc(1)", 400,0,16000,400,0,16000);
  TH2D *c32 = new TH2D("c32", "hgfc(3) vs hgbc(2)", 400,0,16000,400,0,16000);
  TH2D *c33 = new TH2D("c33", "hgfc(3) vs hgbc(3)", 400,0,16000,400,0,16000);

  TH1D *t00 = new TH1D("t00", "dt(ns) = hgft    - hgbt   ", 1000,-5000,5000);
  TH1D *t11 = new TH1D("t11", "dt(ns) = hgft(1) - hgbt(1)", 1000,-5000,5000);
  TH1D *t12 = new TH1D("t12", "dt(ns) = hgft(1) - hgbt(2)", 1000,-5000,5000);
  TH1D *t13 = new TH1D("t13", "dt(ns) = hgft(1) - hgbt(3)", 1000,-5000,5000);
  TH1D *t21 = new TH1D("t21", "dt(ns) = hgft(2) - hgbt(1)", 1000,-5000,5000);
  TH1D *t22 = new TH1D("t22", "dt(ns) = hgft(2) - hgbt(2)", 1000,-5000,5000);
  TH1D *t23 = new TH1D("t23", "dt(ns) = hgft(2) - hgbt(3)", 1000,-5000,5000);
  TH1D *t31 = new TH1D("t31", "dt(ns) = hgft(3) - hgbt(1)", 1000,-5000,5000);
  TH1D *t32 = new TH1D("t32", "dt(ns) = hgft(3) - hgbt(2)", 1000,-5000,5000);
  TH1D *t33 = new TH1D("t33", "dt(ns) = hgft(3) - hgbt(3)", 1000,-5000,5000);

  TH1D *g0  = new TH1D("g0", "size of decay pixels satisfied geo + time", 100,0,100);
  TH2D *g1  = new TH2D("sum", "summary", 4000,0,16000,300,0,300);
  TH2D *g2  = new TH2D("sum_cal", "summary_cal", 4000,0,16000,300,0,300);

  for(int runnum=1040;runnum<1050;runnum++){ 
    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);  

    long nentries = chan->GetEntries();
    long x = 0;


    std::pair<int,int> piximp;
    std::vector<double> decaytime;
    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0){
        piximp = fevent->LGPixel();
      }else{
        int count0 = 0;
        if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
          for(int i=0;i<fevent->HGFSize();i++){
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;
            int    hgfn = fevent->HGF()[i].GetStrip()-1;
            double hgfc = fevent->HGF()[i].GetCharge();
            double hgfe = fevent->HGF()[i].GetEnergy();
            double hgft = fevent->HGF()[i].GetTimestamp();
            if(fevent->Range(hgfn,0,15))  f1 = true;
            if(fevent->Range(hgfn,16,31)) f2 = true;
            if(fevent->Range(hgfn,32,39)) f3 = true;
            g1->Fill(hgfc, hgfn);
            g2->Fill(hgfe, hgfn);
            for(int j=0;j<fevent->HGBSize();j++){
              bool b1 = false;
              bool b2 = false;
              bool b3 = false;
              int    hgbn = fevent->HGB()[j].GetStrip()-1;
              double hgbc = fevent->HGB()[j].GetCharge();
              double hgbe = fevent->HGB()[j].GetEnergy();
              double hgbt = fevent->HGB()[j].GetTimestamp();
              if(fevent->Range(hgbn,0,15))  b1 = true;
              if(fevent->Range(hgbn,16,31)) b2 = true;
              if(fevent->Range(hgbn,32,39)) b3 = true;

              double dt = hgft - hgbt;
              if(dt>=30 && dt<140){ 
                if(fabs(hgfn-piximp.first)<=1 && fabs(hgbn-piximp.second)<=1){
                  count0++;
                  c00->Fill(hgfc,hgbc);
                  t00->Fill((hgft-hgbt)); 
                  if(f1 && b1){ c11->Fill(hgfc, hgbc); t11->Fill((hgft-hgbt)); } 
                  if(f1 && b2){ c12->Fill(hgfc, hgbc); t12->Fill((hgft-hgbt)); } 
                  if(f1 && b3){ c13->Fill(hgfc, hgbc); t13->Fill((hgft-hgbt)); } 
                  if(f2 && b1){ c21->Fill(hgfc, hgbc); t21->Fill((hgft-hgbt)); } 
                  if(f2 && b2){ c22->Fill(hgfc, hgbc); t22->Fill((hgft-hgbt)); } 
                  if(f2 && b3){ c23->Fill(hgfc, hgbc); t23->Fill((hgft-hgbt)); } 
                  if(f3 && b1){ c31->Fill(hgfc, hgbc); t31->Fill((hgft-hgbt)); } 
                  if(f3 && b2){ c32->Fill(hgfc, hgbc); t32->Fill((hgft-hgbt)); } 
                  if(f3 && b3){ c33->Fill(hgfc, hgbc); t33->Fill((hgft-hgbt)); } 
                  g1->Fill(hgbc, hgbn+80);
                  g2->Fill(hgbe, hgbn+80);
                }
              }
            }
          }
        }
        g0->Fill(count0);
      }

      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum, x,nentries);
        fflush(stdout);
      }

    }

    printf("run%i.root: on entry %lu / %lu \n", runnum,x, nentries);
  }
  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2.root","recreate");
  //TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2_blobs.root","recreate");
  c00->Write();
  c11->Write();
  c12->Write();
  c13->Write();
  c21->Write();
  c22->Write();
  c23->Write();
  c31->Write();
  c32->Write();
  c33->Write();
  t00->Write();
  t11->Write();
  t12->Write();
  t13->Write();
  t21->Write();
  t22->Write();
  t23->Write();
  t31->Write();
  t32->Write();
  t33->Write();
  g0 ->Write();
  g1 ->Write();
  g2 ->Write();
  newf->Close();
}







void decaysort_blob(){
  
  TChannel::ReadDetMapFile();
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }
  
    
  TH2D *ct00[6];
  TH2D *ct11[6];
  TH2D *ct12[6];
  TH2D *ct13[6];
  TH2D *ct21[6];
  TH2D *ct22[6];
  TH2D *ct23[6];
  TH2D *ct31[6];
  TH2D *ct32[6];
  TH2D *ct33[6];
  TH2D *c00[6];
  TH2D *c11[6];
  TH2D *c12[6];
  TH2D *c13[6];
  TH2D *c21[6];
  TH2D *c22[6];
  TH2D *c23[6];
  TH2D *c31[6];
  TH2D *c32[6];
  TH2D *c33[6];
  TH1D *g0[6];

  TH2D *h11 = new TH2D("h11","HGFE(1) vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h12 = new TH2D("h12","HGFE(2) vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h13 = new TH2D("h13","HGFE(3) vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h21 = new TH2D("h21","HGBE(1)vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h22 = new TH2D("h22","HGBE(2)vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h23 = new TH2D("h23","HGBE(3)vs dt(hgft - hgbt)"                , 1000,-5000,5000,400,0,16000);
  TH2D *h31 = new TH2D("h31","HGFE(1)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h32 = new TH2D("h32","HGFE(2)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h33 = new TH2D("h33","HGFE(3)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h41 = new TH2D("h41","HGBE(1)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h42 = new TH2D("h42","HGBE(2)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h43 = new TH2D("h43","HGBE(3)vs dt(hgft - hgbt) within dtgate"  , 1000,-5000,5000,400,0,16000);
  TH2D *h51 = new TH2D("h51","HGFE(1)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);
  TH2D *h52 = new TH2D("h52","HGFE(2)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);
  TH2D *h53 = new TH2D("h53","HGFE(3)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);
  TH2D *h61 = new TH2D("h61","HGBE(1)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);
  TH2D *h62 = new TH2D("h62","HGBE(2)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);
  TH2D *h63 = new TH2D("h63","HGBE(3)vs dt(hgft - hgbt) with good pixel", 1000,-5000,5000,400,0,16000);


  for(int i=0;i<6;i++){
    c00[i] = new TH2D(Form("c00_%i",i), Form("hgfc    vs hgbc    in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c11[i] = new TH2D(Form("c11_%i",i), Form("hgfc(1) vs hgbc(1) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c12[i] = new TH2D(Form("c12_%i",i), Form("hgfc(1) vs hgbc(2) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c13[i] = new TH2D(Form("c13_%i",i), Form("hgfc(1) vs hgbc(3) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c21[i] = new TH2D(Form("c21_%i",i), Form("hgfc(2) vs hgbc(1) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c22[i] = new TH2D(Form("c22_%i",i), Form("hgfc(2) vs hgbc(2) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c23[i] = new TH2D(Form("c23_%i",i), Form("hgfc(2) vs hgbc(3) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c31[i] = new TH2D(Form("c31_%i",i), Form("hgfc(3) vs hgbc(1) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c32[i] = new TH2D(Form("c32_%i",i), Form("hgfc(3) vs hgbc(2) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    c33[i] = new TH2D(Form("c33_%i",i), Form("hgfc(3) vs hgbc(3) in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct00[i] = new TH2D(Form("ct00_%i",i), Form("hgfc    vs hgbc    gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct11[i] = new TH2D(Form("ct11_%i",i), Form("hgfc(1) vs hgbc(1) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct12[i] = new TH2D(Form("ct12_%i",i), Form("hgfc(1) vs hgbc(2) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct13[i] = new TH2D(Form("ct13_%i",i), Form("hgfc(1) vs hgbc(3) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct21[i] = new TH2D(Form("ct21_%i",i), Form("hgfc(2) vs hgbc(1) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct22[i] = new TH2D(Form("ct22_%i",i), Form("hgfc(2) vs hgbc(2) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct23[i] = new TH2D(Form("ct23_%i",i), Form("hgfc(2) vs hgbc(3) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct31[i] = new TH2D(Form("ct31_%i",i), Form("hgfc(3) vs hgbc(1) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct32[i] = new TH2D(Form("ct32_%i",i), Form("hgfc(3) vs hgbc(2) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    ct33[i] = new TH2D(Form("ct33_%i",i), Form("hgfc(3) vs hgbc(3) gated time in %s",blob[i]->GetName()), 500,0,25000,500,0,25000);
    g0[i]  = new TH1D(Form("g0_%i",i) , Form("size of decay pixels satisfied geo + time in %s", blob[i]->GetName()), 100,0,100);
  }

  TH2D *g1  = new TH2D("sum", "summary", 4000,0,16000,300,0,300);
  TH2D *g2  = new TH2D("sum_cal", "summary_cal", 4000,0,16000,300,0,300);

  for(int runnum=1040;runnum<1050;runnum++){
    std::vector<double> tofpar= TOFCorrection::Get()->ReadFile(runnum);
    std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    //chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/decay/highgain/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);  

    long nentries = chan->GetEntries();
    long x = 0;

    std::pair<int,int> piximp;
    double tof, pin1e;
    bool flag_blob[6];
    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0){
        for(int i=0;i<6;i++) { flag_blob[i] = false; }
        piximp = fevent->LGPixel();
        pin1e = fevent->Pin1E();
        tof   = fevent->I2S();
        tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
        for(int i=0;i<6;i++){
          if(blob[i]->IsInside(tof,pin1e)){
            flag_blob[i] = true;
          }
        }
      }else{
        int count0 = 0;
        if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
          for(int i=0;i<fevent->HGFSize();i++){
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;
            int    hgfn = fevent->HGF()[i].GetStrip()-1;
            double hgfc = fevent->HGF()[i].GetCharge();
            double hgfe = fevent->HGF()[i].GetEnergy();
            double hgft = fevent->HGF()[i].GetTimestamp();
            if(hgfc<10 || hgfc>25000) continue;
            if(fevent->Range(hgfn,0,15))  f1 = true;
            if(fevent->Range(hgfn,16,31)) f2 = true;
            if(fevent->Range(hgfn,32,39)) f3 = true;
            g1->Fill(hgfc, hgfn);
            g2->Fill(hgfe, hgfn);
            for(int j=0;j<fevent->HGBSize();j++){
              bool b1 = false;
              bool b2 = false;
              bool b3 = false;
              int    hgbn = fevent->HGB()[j].GetStrip()-1;
              double hgbc = fevent->HGB()[j].GetCharge();
              double hgbe = fevent->HGB()[j].GetEnergy();
              double hgbt = fevent->HGB()[j].GetTimestamp();
              if(hgbc<10 || hgbc>25000) continue;
              if(fevent->Range(hgbn,0,15))  b1 = true;
              if(fevent->Range(hgbn,16,31)) b2 = true;
              if(fevent->Range(hgbn,32,39)) b3 = true;

              double dt = hgft - hgbt;
              if(f1) h11->Fill(dt, hgfc);
              if(f2) h12->Fill(dt, hgfc);
              if(f3) h13->Fill(dt, hgfc);
              if(b1) h21->Fill(dt, hgbc);
              if(b2) h22->Fill(dt, hgbc);
              if(b3) h23->Fill(dt, hgbc);
              for(int m=0;m<6;m++){
                if(!flag_blob[m]) continue;
                c00[m]->Fill(hgfc,hgbc);
                if(f1 && b1){ c11[m]->Fill(hgfc, hgbc); } 
                if(f1 && b2){ c12[m]->Fill(hgfc, hgbc); } 
                if(f1 && b3){ c13[m]->Fill(hgfc, hgbc); } 
                if(f2 && b1){ c21[m]->Fill(hgfc, hgbc); } 
                if(f2 && b2){ c22[m]->Fill(hgfc, hgbc); } 
                if(f2 && b3){ c23[m]->Fill(hgfc, hgbc); } 
                if(f3 && b1){ c31[m]->Fill(hgfc, hgbc); } 
                if(f3 && b2){ c32[m]->Fill(hgfc, hgbc); } 
                if(f3 && b3){ c33[m]->Fill(hgfc, hgbc); } 
              }
              if(dt>=30 && dt<140){ 
                if(f1) h31->Fill(dt, hgfc);
                if(f2) h32->Fill(dt, hgfc);
                if(f3) h33->Fill(dt, hgfc);
                if(b1) h41->Fill(dt, hgbc);
                if(b2) h42->Fill(dt, hgbc);
                if(b3) h43->Fill(dt, hgbc);
                if(fabs(hgfn-piximp.first)<=1 && fabs(hgbn-piximp.second)<=1){
                  if(f1) h51->Fill(dt, hgfc);
                  if(f2) h52->Fill(dt, hgfc);
                  if(f3) h53->Fill(dt, hgfc);
                  if(b1) h61->Fill(dt, hgbc);
                  if(b2) h62->Fill(dt, hgbc);
                  if(b3) h63->Fill(dt, hgbc);
                  count0++;
                  for(int m=0;m<6;m++){
                    if(!flag_blob[m]) continue;
                    ct00[m]->Fill(hgfc,hgbc);
                    if(f1 && b1){ ct11[m]->Fill(hgfc, hgbc); } 
                    if(f1 && b2){ ct12[m]->Fill(hgfc, hgbc); } 
                    if(f1 && b3){ ct13[m]->Fill(hgfc, hgbc); } 
                    if(f2 && b1){ ct21[m]->Fill(hgfc, hgbc); } 
                    if(f2 && b2){ ct22[m]->Fill(hgfc, hgbc); } 
                    if(f2 && b3){ ct23[m]->Fill(hgfc, hgbc); } 
                    if(f3 && b1){ ct31[m]->Fill(hgfc, hgbc); } 
                    if(f3 && b2){ ct32[m]->Fill(hgfc, hgbc); } 
                    if(f3 && b3){ ct33[m]->Fill(hgfc, hgbc); } 
                  }
                  g1->Fill(hgbc, hgbn+80);
                  g2->Fill(hgbe, hgbn+80);
                }
              }
            }
          }
        }
        for(int m=0;m<6;m++){
          if(!flag_blob[m]) continue;
          g0[m]->Fill(count0);
        }
      }

      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum, x,nentries);
        fflush(stdout);
      }

    }

    printf("run%i.root: on entry %lu / %lu \n", runnum,x, nentries);
  }
  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2_blobs.root","recreate");
  for(int m=0;m<6;m++){
    c00[m]->Write();
    c11[m]->Write();
    c12[m]->Write();
    c13[m]->Write();
    c21[m]->Write();
    c22[m]->Write();
    c23[m]->Write();
    c31[m]->Write();
    c32[m]->Write();
    c33[m]->Write();
    ct00[m]->Write();
    ct11[m]->Write();
    ct12[m]->Write();
    ct13[m]->Write();
    ct21[m]->Write();
    ct22[m]->Write();
    ct23[m]->Write();
    ct31[m]->Write();
    ct32[m]->Write();
    ct33[m]->Write();
    g0[m] ->Write();
  }
  g1 ->Write();
  g2 ->Write();
  h11 ->Write();
  h12 ->Write();
  h13 ->Write();
  h21 ->Write();
  h22 ->Write();
  h23 ->Write();
  h31 ->Write();
  h32 ->Write();
  h33 ->Write();
  h41 ->Write();
  h42 ->Write();
  h43 ->Write();
  h51 ->Write();
  h52 ->Write();
  h53 ->Write();
  h61 ->Write();
  h62 ->Write();
  h63 ->Write();
  newf->Close();
}





void decay_E_blob(){


  TChannel::ReadDetMapFile();
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }

  TH2D *ef[6];
  TH2D *ef1[6];
  TH2D *ef2[6];
  TH2D *ef3[6];
  TH2D *eb[6];
  TH2D *eb1[6];
  TH2D *eb2[6];
  TH2D *eb3[6];
  TH2D *de[6];
  TH2D *de11[6];
  TH2D *de12[6];
  TH2D *de13[6];
  TH2D *de21[6];
  TH2D *de22[6];
  TH2D *de23[6];
  TH2D *de31[6];
  TH2D *de32[6];
  TH2D *de33[6];
  TH2D *hf[6];
  TH2D *hf1[6];
  TH2D *hf2[6];
  TH2D *hf3[6];
  TH2D *hb[6];
  TH2D *hb1[6];
  TH2D *hb2[6];
  TH2D *hb3[6];

  TH2D *le11[6];
  TH2D *sum[6];

  TH2D *sf[6];
  TH2D *sb[6];
  TH2D *sf0[6];
  TH2D *sb0[6];


  for(int i=0;i<6;i++){
    sum[i] = new TH2D(Form("sum_%i",i), Form("summary charge from decay in %s", blob[i]->GetName()), 4000,0,36000,300,0,300);
    ef[i]  = new TH2D(Form("ef_%i",i) , Form("HGFMAx.E vs dt(hgft-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    eb[i]  = new TH2D(Form("eb_%i",i) , Form("HGBMax.E vs dt(hgft-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    ef1[i] = new TH2D(Form("ef1_%i",i), Form("HGFMax(1).E vs dt(hgft(1)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    ef2[i] = new TH2D(Form("ef2_%i",i), Form("HGFMax(2).E vs dt(hgft(2)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    ef3[i] = new TH2D(Form("ef3_%i",i), Form("HGFMax(3).E vs dt(hgft(3)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    eb1[i] = new TH2D(Form("eb1_%i",i), Form("HGBMax(1).E vs dt(hgft(1)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    eb2[i] = new TH2D(Form("eb2_%i",i), Form("HGBMax(2).E vs dt(hgft(2)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);
    eb3[i] = new TH2D(Form("eb3_%i",i), Form("HGBMax(3).E vs dt(hgft(3)-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,400,0,16000);

    de[i]   = new TH2D(Form("de_%i",i)  ,Form("Dec: dE vs dt(hgft-hgbt) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de11[i] = new TH2D(Form("de11_%i",i),Form("Dec: dE vs dt(hgft(1)-hgbt(1)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de12[i] = new TH2D(Form("de12_%i",i),Form("Dec: dE vs dt(hgft(1)-hgbt(2)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de13[i] = new TH2D(Form("de13_%i",i),Form("Dec: dE vs dt(hgft(1)-hgbt(3)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de21[i] = new TH2D(Form("de21_%i",i),Form("Dec: dE vs dt(hgft(2)-hgbt(1)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de22[i] = new TH2D(Form("de22_%i",i),Form("Dec: dE vs dt(hgft(2)-hgbt(2)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de23[i] = new TH2D(Form("de23_%i",i),Form("Dec: dE vs dt(hgft(2)-hgbt(3)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de31[i] = new TH2D(Form("de31_%i",i),Form("Dec: dE vs dt(hgft(3)-hgbt(1)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de32[i] = new TH2D(Form("de32_%i",i),Form("Dec: dE vs dt(hgft(3)-hgbt(2)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);
    de33[i] = new TH2D(Form("de33_%i",i),Form("Dec: dE vs dt(hgft(3)-hgbt(3)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);

    hf[i]   = new TH2D(Form("hf_%i",i) ,Form("dE vs dt(hgft - sssdt) in %s",blob[i]->GetName())   ,1000,-5000,5000,800,-16000,16000);
    hf1[i]  = new TH2D(Form("hf1_%i",i),Form("dE vs dt(hgft(1) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);
    hf2[i]  = new TH2D(Form("hf2_%i",i),Form("dE vs dt(hgft(2) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);
    hf3[i]  = new TH2D(Form("hf3_%i",i),Form("dE vs dt(hgft(3) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);
    hb[i]   = new TH2D(Form("hb_%i",i) ,Form("dE vs dt(hgbt - sssdt) in %s",blob[i]->GetName())   ,1000,-5000,5000,800,-16000,16000);
    hb1[i]  = new TH2D(Form("hb1_%i",i),Form("dE vs dt(hgbt(1) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);
    hb2[i]  = new TH2D(Form("hb2_%i",i),Form("dE vs dt(hgbt(2) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);
    hb3[i]  = new TH2D(Form("hb3_%i",i),Form("dE vs dt(hgbt(3) - sssdt) in %s",blob[i]->GetName()),1000,-5000,5000,800,-16000,16000);

    le11[i] = new TH2D(Form("le11_%i",i),Form("Imp: dE vs dt(hgft(1)-hgbt(1)) in %s",blob[i]->GetName()), 1000,-5000,5000,800,-16000,16000);

    sf0[i]  = new TH2D(Form("sf0_%i",i), Form("HGFSize vs HGFMax.C %s",blob[i]->GetName()), 40,0,40, 4000,0,36000);
    sb0[i]  = new TH2D(Form("sb0_%i",i), Form("HGBSize vs HGBMax.C %s",blob[i]->GetName()), 40,0,40, 4000,0,36000);
    sf[i]   = new TH2D(Form("sf_%i",i) , Form("HGFSize[10,25000] vs HGFMax.C %s",blob[i]->GetName()), 40,0,40, 4000,0,36000);
    sb[i]   = new TH2D(Form("sb_%i",i) , Form("HGBSize[10,25000] vs HGBMax.C %s",blob[i]->GetName()), 40,0,40, 4000,0,36000);
  }

  
  int count0 = 0;
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  for(int runnum=1040;runnum<1050;runnum++){
    std::vector<double> tofpar= TOFCorrection::Get()->ReadFile(runnum);
    std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);  

    long nentries = chan->GetEntries();
    long x = 0;

    std::pair<int,int> piximp;
    double tof, pin1e;
    bool flag_blob[6];
    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0){
        if(fevent->LGFSize()>10) continue;
        for(int i=0;i<6;i++) { flag_blob[i] = false; }
        piximp = fevent->LGPixel();
        pin1e = fevent->Pin1E();
        tof   = fevent->I2S();
        tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
        for(int i=0;i<6;i++){
          if(blob[i]->IsInside(tof,pin1e)){
            flag_blob[i] = true;
          }
        }
      }else{
        if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
          count0++;
          if(fevent->LGFSize()>0 || fevent->LGBSize()>0) { count1++; continue;}
          for(int m=0;m<6;m++){
            if(!flag_blob[m]) continue;
            sf0[m]->Fill(fevent->HGFSize(),fevent->HGFMax().GetCharge());
            sb0[m]->Fill(fevent->HGBSize(),fevent->HGBMax().GetCharge());
            if(fevent->HGFMax().GetCharge()>=10 && fevent->HGFMax().GetCharge()<25000){
              sf[m] ->Fill(fevent->HGFSize(10,25000), fevent->HGFMax().GetCharge());
            }
            if(fevent->HGBMax().GetCharge()>=10 && fevent->HGBMax().GetCharge()<25000){
              sb[m] ->Fill(fevent->HGBSize(10,25000), fevent->HGBMax().GetCharge());
            }
            if(fevent->HGFSize(10,25000)==0 && fevent->HGBSize(10,25000)==0) count2++;
            if(fevent->HGFSize(10,25000)>0 && fevent->HGBSize(10,25000)>0){
              if(fevent->HGFSize(10,25000)<5 && fevent->HGBSize(10,25000)<5)
                count3++;
            } 
            for(auto it:fevent->fHits){
              sum[m]->Fill(it.GetCharge(), it.GetNumber());
            }
          }
          bool f1, f2, f3, b1, b2, b3;
    
          //double cmax = 0;
          //int hgfn, hgbn;
          //double hgft, hgfc, hgbt, hgbc;
          //for(int i=0;i<fevent->HGFSize();i++){
          //  double tempc = fevent->HGF()[i].GetCharge();
          //  if(tempc<80 || tempc>25000) continue;
          //  if(tempc>cmax){
          //    hgfn = fevent->HGF()[i].GetStrip()-1;
          //    hgft = fevent->HGF()[i].GetTimestamp();
          //    hgfc = fevent->HGF()[i].GetCharge();
          //    cmax = tempc;
          //  } 
          //}
          //if(cmax==0) {count2++; continue;}
          //cmax=0;
          //for(int i=0;i<fevent->HGBSize();i++){
          //  double tempc = fevent->HGB()[i].GetCharge();
          //  if(tempc<80 || tempc>25000) continue;
          //  if(tempc>cmax){
          //    hgbn = fevent->HGB()[i].GetStrip()-1;
          //    hgbt = fevent->HGB()[i].GetTimestamp();
          //    hgbc = fevent->HGB()[i].GetCharge();
          //    cmax = tempc;
          //  } 
          //}
          //if(cmax==0) {count2++; continue;}
          //if(hgbc>3200 && hgbc<5000) continue;
          int    hgfn = fevent->HGFMax().GetStrip()-1;
          double hgft = fevent->HGFMax().GetTimestamp();
          double hgfc = fevent->HGFMax().GetCharge();
          int    hgbn = fevent->HGBMax().GetStrip()-1;
          double hgbt = fevent->HGBMax().GetTimestamp();
          double hgbc = fevent->HGBMax().GetCharge();
          if(fevent->Range(hgfn,0,15))  f1 = true;
          if(fevent->Range(hgfn,16,31)) f2 = true;
          if(fevent->Range(hgfn,32,39)) f3 = true;
          if(fevent->Range(hgbn,0,15))  b1 = true;
          if(fevent->Range(hgbn,16,31)) b2 = true;
          if(fevent->Range(hgbn,32,39)) b3 = true;
          double dt = hgft - hgbt;
          double dc = hgfc - hgbc;
          for(int m=0;m<6;m++){
            if(!flag_blob[m]) continue;
            ef[m]->Fill(dt, hgfc);    
            eb[m]->Fill(dt, hgbc);    
            if(f1) { ef1[m]->Fill(dt,hgfc); }
            if(f2) { ef2[m]->Fill(dt,hgfc); }
            if(f3) { ef3[m]->Fill(dt,hgfc); }
            if(b1) { eb1[m]->Fill(dt,hgbc); }
            if(b2) { eb2[m]->Fill(dt,hgbc); }
            if(b3) { eb3[m]->Fill(dt,hgbc); }

            de[m]->Fill(dt,dc);
            if(f1 && b1) { de11[m]->Fill(dt,dc); }    
            if(f1 && b2) { de12[m]->Fill(dt,dc); }    
            if(f1 && b3) { de13[m]->Fill(dt,dc); }    
            if(f2 && b1) { de21[m]->Fill(dt,dc); }    
            if(f2 && b2) { de22[m]->Fill(dt,dc); }    
            if(f2 && b3) { de23[m]->Fill(dt,dc); }    
            if(f3 && b1) { de31[m]->Fill(dt,dc); }    
            if(f3 && b2) { de32[m]->Fill(dt,dc); }    
            if(f3 && b3) { de33[m]->Fill(dt,dc); }    
          }    

          //for(int i=0;i<fevent->HGFSize();i++){
          //  bool f1 = false;
          //  bool f2 = false;
          //  bool f3 = false;
          //  int    hgfn = fevent->HGF()[i].GetStrip()-1;
          //  double hgfc = fevent->HGF()[i].GetCharge();
          //  double hgfe = fevent->HGF()[i].GetEnergy();
          //  double hgft = fevent->HGF()[i].GetTimestamp();
          //  if(hgfc<80 || hgfc>25000) continue;
          //  if(fevent->Range(hgfn,0,15))  f1 = true;
          //  if(fevent->Range(hgfn,16,31)) f2 = true;
          //  if(fevent->Range(hgfn,32,39)) f3 = true;
          //  for(int j=0;j<fevent->HGBSize();j++){
          //    bool b1 = false;
          //    bool b2 = false;
          //    bool b3 = false;
          //    int    hgbn = fevent->HGB()[j].GetStrip()-1;
          //    double hgbc = fevent->HGB()[j].GetCharge();
          //    double hgbe = fevent->HGB()[j].GetEnergy();
          //    double hgbt = fevent->HGB()[j].GetTimestamp();
          //    if(hgbc<80 || hgbc>25000) continue;
          //    if(fevent->Range(hgbn,0,15))  b1 = true;
          //    if(fevent->Range(hgbn,16,31)) b2 = true;
          //    if(fevent->Range(hgbn,32,39)) b3 = true;
          //    double dt = hgft - hgbt;
          //    double dc = hgfc - hgbc;    
          //    for(int m=0;m<6;m++){
          //      if(!flag_blob[m]) continue;
          //      ef[m]->Fill(dt, hgfc);
          //      eb[m]->Fill(dt, hgbc);
          //      if(f1) {ef1[m]->Fill(dt, hgfc); }
          //      if(f2) {ef2[m]->Fill(dt, hgfc); }
          //      if(f3) {ef3[m]->Fill(dt, hgfc); }
          //      if(b1) {eb1[m]->Fill(dt, hgbc); }
          //      if(b2) {eb2[m]->Fill(dt, hgbc); }
          //      if(b3) {eb3[m]->Fill(dt, hgbc); }
          //      
          //      de[m]->Fill(dt, dc);
          //      if(f1 && b1) { de11[m]->Fill(dt,dc); }
          //      if(f1 && b2) { de12[m]->Fill(dt,dc); }
          //      if(f1 && b3) { de13[m]->Fill(dt,dc); }
          //      if(f2 && b1) { de21[m]->Fill(dt,dc); }
          //      if(f2 && b2) { de22[m]->Fill(dt,dc); }
          //      if(f2 && b3) { de23[m]->Fill(dt,dc); }
          //      if(f3 && b1) { de31[m]->Fill(dt,dc); }
          //      if(f3 && b2) { de32[m]->Fill(dt,dc); }
          //      if(f3 && b3) { de33[m]->Fill(dt,dc); }
          //    }
          //  }// HGB loop over
          //}// HGF loop over
        
          //for(int i=0;i<fevent->SSSDSize();i++){
          //  double sssdc = fevent->SSSD()[i].GetCharge();
          //  double sssdt = fevent->SSSD()[i].GetTimestamp();
          //  if(sssdc<80 || sssdc>25000) continue;
          //  for(int j=0;j<fevent->HGFSize();j++){
          //    bool f1 = false;
          //    bool f2 = false;
          //    bool f3 = false;
          //    int    hgfn = fevent->HGF()[j].GetStrip()-1;
          //    double hgfc = fevent->HGF()[j].GetCharge();
          //    double hgft = fevent->HGF()[j].GetTimestamp();
          //    if(hgfc<80 || hgfc>25000) continue;
          //    if(fevent->Range(hgfn,0,15))  f1 = true;
          //    if(fevent->Range(hgfn,16,31)) f2 = true;
          //    if(fevent->Range(hgfn,32,39)) f3 = true;
          //    double dt = hgft - sssdt;
          //    double dc = hgfc - sssdc;
          //    for(int m=0;m<6;m++){
          //      if(!flag_blob[m]) continue;
          //      hf[m]->Fill(dt,dc);
          //      if(f1) { hf1[m]->Fill(dt,dc); }
          //      if(f2) { hf2[m]->Fill(dt,dc); }
          //      if(f3) { hf3[m]->Fill(dt,dc); }
          //    }
          //  }  
          //  for(int j=0;j<fevent->HGBSize();j++){
          //    bool b1 = false;
          //    bool b2 = false;
          //    bool b3 = false;
          //    int    hgbn = fevent->HGB()[j].GetStrip()-1;
          //    double hgbc = fevent->HGB()[j].GetCharge();
          //    double hgbt = fevent->HGB()[j].GetTimestamp();
          //    if(hgbc<80 || hgbc>25000) continue;
          //    if(fevent->Range(hgbn,0,15))  b1 = true;
          //    if(fevent->Range(hgbn,16,31)) b2 = true;
          //    if(fevent->Range(hgbn,32,39)) b3 = true;
          //    double dt = hgbt - sssdt;
          //    double dc = hgbc - sssdc;
          //    for(int m=0;m<6;m++){
          //      if(!flag_blob[m]) continue;
          //      hb[m]->Fill(dt,dc);
          //      if(b1) { hb1[m]->Fill(dt,dc); }
          //      if(b2) { hb2[m]->Fill(dt,dc); }
          //      if(b3) { hb3[m]->Fill(dt,dc); }
          //    }
          //  }  
          //}// SSSD loop over 

        }
      }

      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum, x,nentries);
        fflush(stdout);
      }

    }

    printf("run%i.root: on entry %lu / %lu \n", runnum,x, nentries);
  }

  std::cout<<"(dec) hgf && hgb = " << count0 << std::endl;
  std::cout<<"(dec) hgf && hgb && (lgf || lgb) = " << count1 << std::endl;
  std::cout<<"(dec) hgf && hgb && (lgf || lgb) && no HG.C in [10,25000] = " << count2 << std::endl;
  std::cout<<"(dec) hgf && hgb && (no lgf && no lgb) && no HGF&B.size(10,25000) = [1,2] = " << count3 << std::endl;


  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2_decay_nolg_maxC_blobs.root","recreate");
  for(int i=0;i<6;i++){
    sum[i] ->Write();
    ef[i]  ->Write();
    eb[i]  ->Write();
    ef1[i] ->Write();
    ef2[i] ->Write();
    ef3[i] ->Write();
    eb1[i] ->Write();
    eb2[i] ->Write();
    eb3[i] ->Write();

    de[i]  ->Write();
    de11[i]->Write();
    de12[i]->Write();
    de13[i]->Write();
    de21[i]->Write();
    de22[i]->Write();
    de23[i]->Write();
    de31[i]->Write();
    de32[i]->Write();
    de33[i]->Write();

    sf0[i] ->Write();
    sb0[i] ->Write();
    sf[i]  ->Write();
    sb[i]  ->Write();
    //hf[i]->Write();
    //hf1[i]->Write();
    //hf2[i]->Write();
    //hf3[i]->Write();
    //hb[i]->Write();
    //hb1[i]->Write();
    //hb2[i]->Write();
    //hb3[i]->Write();

    //le11[i] ->Write();
  }
  newf->Close();

}


void decay_Tfluctuation_check(){

  TChannel::ReadDetMapFile();
  
  TH2D *hf1  = new TH2D("hf1"  , "HGF(1) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hf2  = new TH2D("hf2"  , "HGF(2) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hf3  = new TH2D("hf3"  , "HGF(3) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hb1  = new TH2D("hb1"  , "HGB(1) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hb2  = new TH2D("hb2"  , "HGB(2) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hb3  = new TH2D("hb3"  , "HGB(3) Charge vs running time", 2000,0,2000000,500,0,25000);
  TH2D *hsssd= new TH2D("hsssd", "SSSD Charge vs running time"  , 2000,0,2000000,500,0,25000);
  TH2D *sum = new TH2D("sum", "summary charge from decay", 4000,0,36000,300,0,300);

  long n0 = 0;

  for(int runnum=1040;runnum<1050;runnum++){ 
    TChain *chan = new TChain("event"); 
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);  

    long nentries = chan->GetEntries();
    long x = 0;

    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0) continue;
      for(auto it:fevent->fHits){
        sum->Fill(it.GetCharge(), it.GetNumber());
        bool f1, f2, f3, b1, b2, b3,sssd;
        if(it.GetCharge()>25000) continue;
        if(fevent->Range(it.GetNumber(),0,15))    f1   = true;
        if(fevent->Range(it.GetNumber(),16,31))   f2   = true;
        if(fevent->Range(it.GetNumber(),32,39))   f3   = true;
        if(fevent->Range(it.GetNumber(),80,95))   b1   = true;
        if(fevent->Range(it.GetNumber(),96,111))  b2   = true;
        if(fevent->Range(it.GetNumber(),112,119)) b3   = true;
        if(fevent->Range(it.GetNumber(),160,175)) sssd = true;
        if(f1)  { hf1  ->Fill(x+n0, it.GetCharge()); }
        if(f2)  { hf2  ->Fill(x+n0, it.GetCharge()); }
        if(f3)  { hf3  ->Fill(x+n0, it.GetCharge()); }
        if(b1)  { hb1  ->Fill(x+n0, it.GetCharge()); }
        if(b2)  { hb2  ->Fill(x+n0, it.GetCharge()); }
        if(b3)  { hb3  ->Fill(x+n0, it.GetCharge()); }
        if(sssd){ hsssd->Fill(x+n0, it.GetCharge()); }
      }// fHits loop end
      
      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r", runnum,x,nentries);
        fflush(stdout);
      }
    }
    n0 += nentries;
    printf("run%i.root: on entry %lu / %lu \n", runnum,x,nentries);
  }

  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/decay/correlation2_decay_tfluctuationcheck.root","recreate");
  hf1  ->Write();
  hf2  ->Write();
  hf3  ->Write();
  hb1  ->Write();
  hb2  ->Write();
  hb3  ->Write();
  hsssd->Write();
  
  newf->Close();
}



void decay_entry_check(){


  TChannel::ReadDetMapFile();
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }
  
  int runnum = 1040;
  std::vector<double> tofpar= TOFCorrection::Get()->ReadFile(runnum);
  std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);

  long x = 0;
  long nentries = chan->GetEntries();


  double pin1e, pin1t, tof;
  bool flag_blob[6];
  for(x=2e4;x<2.1e4;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->Pin1E()>0){
      if(fevent->LGFSize()>10) continue;
      for(int i=0;i<6;i++) { flag_blob[i] = false; }
      pin1e = fevent->Pin1E();
      pin1t = fevent->Pin1T();
      tof   = fevent->I2S();
      tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
      for(int i=0;i<6;i++){
        if(blob[i]->IsInside(tof,pin1e)){
          flag_blob[i] = true;
        }
      }
    }else{
      if(fevent->LGFSize()>0 || fevent->LGBSize()>0) continue;
      if(fevent->HGFSize(10,25000)<3 || fevent->HGFSize(10,25000)>5) continue;
      if(fevent->HGBSize(10,25000)<3 || fevent->HGBSize(10,25000)>5) continue;
      //if(fevent->HGFMax().GetCharge()<7000 || fevent->HGBMax().GetCharge()<7000) continue;
      bool flag = false;
      for(int m=0;m<6;m++){
        if(flag_blob[m]){
          std::cout<<std::endl;
          std::cout<<blob[m]->GetName()<<"\t"<<"pin1t="<<pin1t<<"\t"<<"dec entry="<<x<<std::endl;
          flag = true;
        }
      }
      if(!flag) continue;
      for(int i=0;i<fevent->HGFSize();i++){
        int    hgfn = fevent->HGF()[i].GetStrip()-1;
        double hgfc = fevent->HGF()[i].GetCharge();
        double hgft = fevent->HGF()[i].GetTimestamp();
        std::cout<<"HGF"<<"\t"<<"n="<<hgfn<<"\t"<<"c="<<hgfc<<"\t"<<"t="<<hgft<<std::endl;
      }
      for(int i=0;i<fevent->HGBSize();i++){
        int    hgbn = fevent->HGB()[i].GetStrip()-1;
        double hgbc = fevent->HGB()[i].GetCharge();
        double hgbt = fevent->HGB()[i].GetTimestamp();
        std::cout<<"HGB"<<"\t"<<"n="<<hgbn<<"\t"<<"c="<<hgbc<<"\t"<<"t="<<hgbt<<std::endl;
      }
    }    
  }  



}




void quick(){


  TChannel::ReadDetMapFile();
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }


  TH1D *t  = new TH1D("t" , "dt(ns) = hgft -hgbt"                    ,1000,-5000,5000);
  TH1D *t0 = new TH1D("t0", "dt(ns) = hgft -hgbt in charge[10,25000]",1000,-5000,5000);
  TH1D *t1 = new TH1D("t1", "dt(ns) = hgft -hgbt with maxC"          ,1000,-5000,5000);
  
  TH2D *h  = new TH2D("pid" , "pid"              , 800,10000,18000,300,4200,7200);
  TH2D *h1 = new TH2D("pid1", "pid in delay gate", 800,10000,18000,300,4200,7200);

  for(int runnum=1040;runnum<1050;runnum++){
    std::vector<double> tofpar= TOFCorrection::Get()->ReadFile(runnum);
    std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/correlation2/corlgeo/event%i*.root",runnum));
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);  

    long nentries = chan->GetEntries();
    long x = 0;

    
    double tof, pin1e, pin1t;
    bool flag_blob[6];
    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      if(fevent->Pin1E()>0){
        if(fevent->LGFSize()>10) continue;
        for(int i=0;i<6;i++) { flag_blob[i] = false; }
        pin1e = fevent->Pin1E();
        pin1t = fevent->Pin1T();
        tof   = fevent->I2S();
        tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0]+tofpar1[0];
        for(int i=0;i<6;i++){
          if(blob[i]->IsInside(tof,pin1e)){
            flag_blob[i] = true;
          }
        }
        h->Fill(tof,pin1e);
      }else{
        if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
          if(fevent->LGFSize()>0 || fevent->LGBSize()>0) continue;
          if(fevent->HGFMax().GetEnergy()>0 && fevent->HGBMax().GetEnergy()>0){
            double dt = fevent->HGFMax().GetTimestamp() - fevent->HGBMax().GetTimestamp();
            t1->Fill(dt);
            if(dt>=-220 && dt<-110) h1->Fill(tof, pin1e);
            if(dt>=290  && dt<400) h1->Fill(tof, pin1e);
          }
          for(int i=0;i<fevent->HGFSize();i++){
            int    hgfn = fevent->HGF()[i].GetStrip()-1;
            double hgfc = fevent->HGF()[i].GetCharge();
            double hgft = fevent->HGF()[i].GetTimestamp();
            for(int j=0;j<fevent->HGBSize();j++){
              int    hgbn = fevent->HGB()[j].GetStrip()-1;
              double hgbc = fevent->HGB()[j].GetCharge();
              double hgbt = fevent->HGB()[j].GetTimestamp();
              double dt = hgft - hgbt;
              t->Fill(dt);
              if(hgfc>10 && hgfc<25000 && hgbc>10 && hgbc<25000){
                t0->Fill(dt);
              }
            }
          }
        }
      }
      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum, x,nentries);
        fflush(stdout);
      }
    }
    printf("run%i.root: on entry %lu / %lu \n",runnum, x,nentries);
  }







}



























