

void hist(int runnum = 1038){

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/vetogate.root");
  TCutG *banana1 = (TCutG *)cutf2->Get("banana1");
  TCutG *banana2 = (TCutG *)cutf2->Get("banana2");
 
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");

  
  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  TH2D *h1 = new TH2D("pid1","pid",1000,8000,18000,300,4200,7200);  
  TH2D *h2 = new TH2D("pid2","pid with lgf",1000,8000,18000,300,4200,7200);  
  TH2D *h3 = new TH2D("pid3","pid with lgf+lgb",1000,8000,18000,300,4200,7200);  
  
  TH2D *g1  = new TH2D("hitpad","hitpad on DSSD LG",40,0,40,40,0,40);
  TH2D *g2  = new TH2D("e", "dE(LGFMax) vs E(SSSD)",400,0,36000,400,0,4000);
  TH2D *g21 = new TH2D("e1", "dE(LGFMax1) vs E(SSSD)",400,0,36000,400,0,4000);
  TH2D *g22 = new TH2D("e2", "dE(LGFMax2) vs E(SSSD)",400,0,36000,400,0,4000);
  TH2D *g23 = new TH2D("e3", "dE(LGFMax3) vs E(SSSD)",400,0,36000,400,0,4000);
  TH2D *g3  = new TH2D("elg", "dE(LGFMax) vs E(SSSD_LG)",400,0,36000,400,0,4000);
  TH2D *g31 = new TH2D("elg1", "dE(LGFMax1) vs E(SSSD_LG)",400,0,36000,400,0,4000);
  TH2D *g32 = new TH2D("elg2", "dE(LGFMax2) vs E(SSSD_LG)",400,0,36000,400,0,4000);
  TH2D *g33 = new TH2D("elg3", "dE(LGFMax3) vs E(SSSD_LG)",400,0,36000,400,0,4000);
  TH2D *g4  = new TH2D("elgall", "dE(LGF) vs E(SSSD_LG)",400,0,36000,400,0,4000);

  double time = -1;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    h1->Fill(tof, pin1e);
    if(fevent->LGFSize()==0) continue;
    //printf("entry = %lu\n",x);
    h2->Fill(tof, pin1e);
    if(fevent->LGBSize()>0){
      h3->Fill(tof,pin1e);
      for(int i=0;i<fevent->LGFSize();i++){
        for(int j=0;j<fevent->LGBSize();j++){
          g1->Fill(fevent->LGF().at(i).GetStrip()-1, fevent->LGB().at(j).GetStrip()-1);
        }
      }
    }
    if(fevent->SSSDSize()>0){
      for(int i=0;i<fevent->SSSDSize();i++){
        double sssd  = fevent->SSSD().at(i).GetCharge();
        double sssdn = fevent->SSSD().at(i).GetStrip()-1;
        if(!cutstrip->IsInside(sssdn,fevent->LGFMax().GetStrip()-1)) continue;
        if(fevent->LGFMax1().GetEnergy()>0) g21->Fill(sssd,fevent->LGFMax1().GetEnergy());
        if(fevent->LGFMax2().GetEnergy()>0) g22->Fill(sssd,fevent->LGFMax2().GetEnergy());
        if(fevent->LGFMax3().GetEnergy()>0) g23->Fill(sssd,fevent->LGFMax3().GetEnergy());
        g2->Fill(sssd,fevent->LGFMax().GetEnergy());
      }
    }
    if(fevent->SSSDLGSize()>0){
      for(int i=0;i<fevent->SSSDLGSize();i++){
        double sssdlg = fevent->SSSDLG().at(i).GetCharge();
        if(fevent->LGFMax1().GetEnergy()>0) g31->Fill(sssdlg,fevent->LGFMax1().GetEnergy());
        if(fevent->LGFMax2().GetEnergy()>0) g32->Fill(sssdlg,fevent->LGFMax2().GetEnergy());
        if(fevent->LGFMax3().GetEnergy()>0) g33->Fill(sssdlg,fevent->LGFMax3().GetEnergy());
        g3->Fill(sssdlg,fevent->LGFMax().GetEnergy());
        for(int j=0;j<fevent->LGFSize();j++){
          g4->Fill(sssdlg,fevent->LGF().at(j).GetEnergy());
        }
      }
    }
     

    if((x%5000)==0){
    printf("on entry %lu / %lu \r",x,nentries);
    fflush(stdout);
    } 
  }
  printf("on entry %lu / %lu \n",x,nentries);
  

  TFile *outfile = new TFile(Form("histogram%i.root",runnum),"recreate");
  h1->Write();  
  h2->Write();  
  h3->Write();  
  g1->Write();  
  g21->Write();  
  g22->Write();  
  g23->Write();  
  g31->Write();  
  g32->Write();  
  g33->Write();  
  g2 ->Write();  
  g3 ->Write();  
  g4 ->Write();  
  outfile->Close();

}

void beyondveto(){

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate.root");
  TCutG *ban11 = (TCutG *)cutf2->Get("ban11");
  TCutG *ban12 = (TCutG *)cutf2->Get("ban12");
  TCutG *ban2  = (TCutG *)cutf2->Get("ban2");
  TCutG *ban31 = (TCutG *)cutf2->Get("ban31");
  TCutG *ban32 = (TCutG *)cutf2->Get("ban32");
 
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");

  int runnum[9] = {1030,1031,1032,1033,1034,1035,1036,1037,1039};
 
  TChain *chan = new TChain("event");
  for(int i=0;i<9;i++){
    chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i*.root",runnum[i]));
  }
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  TH2D *h1 = new TH2D("pid1","pid with lgf",1000,8000,18000,300,4200,7200);  
  TH2D *h2 = new TH2D("pid2","pid beyond veto",1000,8000,18000,300,4200,7200);  

  TH2D *g0 = new TH2D("dt0", "dt(ns) LGF - SSSD",1000,-5000,5000,400,0,4000);
  TH2D *g1 = new TH2D("dt1", "dt(ns) PIN1(E) - LGF",1000,-5000,5000,400,2000,10000);
  TH2D *g2 = new TH2D("dt2", "dt(ns) PIN2(E) - LGF",1000,-5000,5000,400,10000,18000);

  TH1D *t0 = new TH1D("DSSDT0", "dt(ns) LGB - LGF", 1000,-5000,5000);
  TH1D *t1 = new TH1D("DSSDT1", "dt(ns) LGBMax - LGFMax", 1000,-5000,5000);
  TH2D *t2 = new TH2D("DSSDT2", "dt(ns) LGB1 - LGB2", 1000,-5000,5000,400,0,4000);
  TH2D *t3 = new TH2D("DSSDT3", "dt(ns) LGF1 - LGF2", 1000,-5000,5000,400,0,4000);

  TH2D *e   = new TH2D("e",  "dE(PIN1) vs E(LGFMax)", 400,0,4000,400,2000,10000);
  TH2D *eb  = new TH2D("eb",  "dE(PIN1) vs E(LGFMax) beyond veto", 400,0,4000,400,2000,10000);
  TH2D *e0  = new TH2D("e0", "dE(PIN2) vs E(LGFMax)", 400,0,4000,400,10000,18000);
  TH2D *e0b = new TH2D("e0b", "dE(PIN2) vs E(LGFMax) beyond veto", 400,0,4000,400,10000,18000);
  TH2D *e11 = new TH2D("e11","LGF1E vs LGB1E", 400,0,4000,400,0,4000);
  TH2D *e12 = new TH2D("e12","LGF1E vs LGB2E", 400,0,4000,400,0,4000);
  TH2D *e13 = new TH2D("e13","LGF1E vs LGB3E", 400,0,4000,400,0,4000);
  TH2D *e21 = new TH2D("e21","LGF2E vs LGB1E", 400,0,4000,400,0,4000);
  TH2D *e22 = new TH2D("e22","LGF2E vs LGB2E", 400,0,4000,400,0,4000);
  TH2D *e23 = new TH2D("e23","LGF2E vs LGB3E", 400,0,4000,400,0,4000);
  TH2D *e31 = new TH2D("e31","LGF3E vs LGB1E", 400,0,4000,400,0,4000);
  TH2D *e32 = new TH2D("e32","LGF3E vs LGB2E", 400,0,4000,400,0,4000);
  TH2D *e33 = new TH2D("e33","LGF3E vs LGB3E", 400,0,4000,400,0,4000);
  TH2D *emax11 = new TH2D("emax11","LGFMax1E vs LGBMax1E", 400,0,4000,400,0,4000);
  TH2D *emax12 = new TH2D("emax12","LGFMax1E vs LGBMax2E", 400,0,4000,400,0,4000);
  TH2D *emax13 = new TH2D("emax13","LGFMax1E vs LGBMax3E", 400,0,4000,400,0,4000);
  TH2D *emax21 = new TH2D("emax21","LGFMax2E vs LGBMax1E", 400,0,4000,400,0,4000);
  TH2D *emax22 = new TH2D("emax22","LGFMax2E vs LGBMax2E", 400,0,4000,400,0,4000);
  TH2D *emax23 = new TH2D("emax23","LGFMax2E vs LGBMax3E", 400,0,4000,400,0,4000);
  TH2D *emax31 = new TH2D("emax31","LGFMax3E vs LGBMax1E", 400,0,4000,400,0,4000);
  TH2D *emax32 = new TH2D("emax32","LGFMax3E vs LGBMax2E", 400,0,4000,400,0,4000);
  TH2D *emax33 = new TH2D("emax33","LGFMax3E vs LGBMax3E", 400,0,4000,400,0,4000);

  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->LGFSize()==0) continue;
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    double lgfmax = fevent->LGFMax().GetEnergy();
    double dt = 0;
    h1->Fill(tof, pin1e);
    e->Fill(lgfmax,pin1e);
    if(fevent->Pin2E()>0){
      e0->Fill(lgfmax, fevent->Pin2E());
    }
    // check veto
    bool veto = false;
    if(fevent->SSSDSize()>0){
      for(int i=0;i<fevent->SSSDSize();i++){
        double sssdn = fevent->SSSD().at(i).GetStrip()-1;
        double sssdc = fevent->SSSD().at(i).GetCharge();
        double sssdt = fevent->SSSD().at(i).GetTimestamp();
        if(cutstrip->IsInside(sssdn, fevent->LGFMax().GetStrip()-1)){
          if(fevent->LGFMax1().GetEnergy()>0){
            if(ban11->IsInside(sssdc, lgfmax) || ban12->IsInside(sssdc, lgfmax)){
              veto = true;
              dt = sssdt - fevent->LGFMax().GetTimestamp();
              g0->Fill(dt,lgfmax); 
            }
          } 
          if(fevent->LGFMax2().GetEnergy()>0){
            if(ban2->IsInside(sssdc, lgfmax)){
              veto = true;
              dt = sssdt - fevent->LGFMax().GetTimestamp();
              g0->Fill(dt,lgfmax); 
            }
          } 
          if(fevent->LGFMax3().GetEnergy()>0){
            if(ban31->IsInside(sssdc, lgfmax) || ban32->IsInside(sssdc, lgfmax)){
              veto = true;
              dt = sssdt - fevent->LGFMax().GetTimestamp();
              g0->Fill(dt,lgfmax); 
            }
          } 
        }
      }
    }// end checking veto

    if(veto) continue;
    eb->Fill(lgfmax,pin1e);
    h2->Fill(tof, pin1e);
    dt = fevent->LGFMax().GetTimestamp() - fevent->Pin1T();
    g1->Fill(dt, fevent->Pin1E());
    if(fevent->Pin2E()>0){
      dt = fevent->LGFMax().GetTimestamp() - fevent->Pin2T();
      g2->Fill(dt, fevent->Pin2E());
      e0b->Fill(lgfmax, fevent->Pin2E());
    }
    if(fevent->LGFSize()>1){
      dt = fevent->LGF().at(0).GetTimestamp() - fevent->LGF().at(1).GetTimestamp();
      t3->Fill(dt,lgfmax);
    }
    if(fevent->LGBSize()>0){
      double lgbmax = fevent->LGBMax().GetEnergy();
      dt = fevent->LGBMax().GetTimestamp() - fevent->LGFMax().GetTimestamp();
      t1->Fill(dt);
      if(fevent->LGBSize()>1){
        dt = fevent->LGB().at(0).GetTimestamp() - fevent->LGB().at(1).GetTimestamp();
        t2->Fill(dt,lgbmax);
      }
      if(fevent->LGFMax1().GetEnergy()>0 && fevent->LGBMax1().GetEnergy()>0) emax11->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax1().GetEnergy()>0 && fevent->LGBMax2().GetEnergy()>0) emax12->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax1().GetEnergy()>0 && fevent->LGBMax3().GetEnergy()>0) emax13->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax2().GetEnergy()>0 && fevent->LGBMax1().GetEnergy()>0) emax21->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax2().GetEnergy()>0 && fevent->LGBMax2().GetEnergy()>0) emax22->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax2().GetEnergy()>0 && fevent->LGBMax3().GetEnergy()>0) emax23->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax3().GetEnergy()>0 && fevent->LGBMax1().GetEnergy()>0) emax31->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax3().GetEnergy()>0 && fevent->LGBMax2().GetEnergy()>0) emax32->Fill(lgfmax,lgbmax);
      if(fevent->LGFMax3().GetEnergy()>0 && fevent->LGBMax3().GetEnergy()>0) emax33->Fill(lgfmax,lgbmax);
      for(int i=0;i<fevent->LGFSize();i++){
        double lgft= fevent->LGF().at(i).GetTimestamp();
        double lgf = fevent->LGF().at(i).GetEnergy();
        int lgfn   = fevent->LGF().at(i).GetNumber();
        for(int j=0;j<fevent->LGBSize();j++){
          double lgbt= fevent->LGB().at(j).GetTimestamp();
          double lgb = fevent->LGB().at(j).GetEnergy();
          int lgbn   = fevent->LGB().at(j).GetNumber();
          if((lgfn>=64 && lgfn<=79) && (lgbn>=144 && lgbn<=159)) e11->Fill(lgf,lgb);
          if((lgfn>=64 && lgfn<=79) && (lgbn>=128 && lgbn<=143)) e12->Fill(lgf,lgb);
          if((lgfn>=64 && lgfn<=79) && (lgbn>=120 && lgbn<=127)) e13->Fill(lgf,lgb);
          if((lgfn>=48 && lgfn<=63) && (lgbn>=144 && lgbn<=159)) e21->Fill(lgf,lgb);
          if((lgfn>=48 && lgfn<=63) && (lgbn>=128 && lgbn<=143)) e22->Fill(lgf,lgb);
          if((lgfn>=48 && lgfn<=63) && (lgbn>=120 && lgbn<=127)) e23->Fill(lgf,lgb);
          if((lgfn>=40 && lgfn<=47) && (lgbn>=144 && lgbn<=159)) e31->Fill(lgf,lgb);
          if((lgfn>=40 && lgfn<=47) && (lgbn>=128 && lgbn<=143)) e32->Fill(lgf,lgb);
          if((lgfn>=40 && lgfn<=47) && (lgbn>=120 && lgbn<=127)) e33->Fill(lgf,lgb);
          dt = lgbt - lgft;
          t0->Fill(dt);
        }
      }
    }

    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
    
  }
  printf("on entry %lu / %lu \n",x,nentries);


  //TFile *outfile = new TFile(Form("histogram%i.root",runnum),"recreate");
  TFile *outfile = new TFile("beyondveto103.root","recreate");
  h1 ->Write();
  h2 ->Write();
  g0 ->Write();
  g1 ->Write();
  g2 ->Write();
  t0 ->Write();
  t1 ->Write();
  t2 ->Write();
  t3 ->Write();   
  e0 ->Write();
  e0b->Write();
  e  ->Write();
  eb ->Write();
  e11->Write();
  e12->Write();
  e13->Write();
  e21->Write();
  e22->Write();
  e23->Write();
  e31->Write();
  e32->Write();
  e33->Write();
  emax11->Write();
  emax12->Write();
  emax13->Write();
  emax21->Write();
  emax22->Write();
  emax23->Write();
  emax31->Write();
  emax32->Write();
  emax33->Write();
  outfile->Close();


}


void beyondveto1(){

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate.root");
  TCutG *ban11 = (TCutG *)cutf2->Get("ban11");
  TCutG *ban12 = (TCutG *)cutf2->Get("ban12");
  TCutG *ban2  = (TCutG *)cutf2->Get("ban2");
  TCutG *ban31 = (TCutG *)cutf2->Get("ban31");
  TCutG *ban32 = (TCutG *)cutf2->Get("ban32");
 
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");
  
  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/DSSDLGE.root");
  TCutG *lge[3];
  for(int i=0;i<3;i++){
    lge[i] = (TCutG *)cutf4->Get(Form("dssdlgfe%i",i+1));
  }


  int runnum[9] = {1030,1031,1032,1033,1034,1035,1036,1037,1039}; 
  TChain *chan = new TChain("event");
  for(int i=0;i<9;i++){
    chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i*.root",runnum[i]));
  }
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  TH2D *h1 = new TH2D("pid1","pid outside veto",1000,8000,18000,300,4200,7200);  
  TH2D *h2 = new TH2D("pid2","pid outside veto inside Egate",1000,8000,18000,300,4200,7200);  

  TH1D *t1 = new TH1D("DSSDT1", "dt(ns) LGB - LGF inside Egate", 1000,-5000,5000);

  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->LGFSize()==0) continue;
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    double lgfmax = fevent->LGFMax().GetEnergy();
    double dt = 0;
    // check veto
    bool veto = false;
    if(fevent->SSSDSize()>0){
      for(int i=0;i<fevent->SSSDSize();i++){
        double sssdn = fevent->SSSD().at(i).GetStrip()-1;
        double sssdc = fevent->SSSD().at(i).GetCharge();
        double sssdt = fevent->SSSD().at(i).GetTimestamp();
        if(cutstrip->IsInside(sssdn, fevent->LGFMax().GetStrip()-1)){
          if(fevent->LGFMax1().GetEnergy()>0){
            if(ban11->IsInside(sssdc, lgfmax) || ban12->IsInside(sssdc, lgfmax)){
              veto = true;
              break;
            }
          } 
          if(fevent->LGFMax2().GetEnergy()>0){
            if(ban2->IsInside(sssdc, lgfmax)){
              veto = true;
              break;
            }
          } 
          if(fevent->LGFMax3().GetEnergy()>0){
            if(ban31->IsInside(sssdc, lgfmax) || ban32->IsInside(sssdc, lgfmax)){
              veto = true;
              break;
            }
          } 
        }
      }
    }// end checking veto


    if(veto) continue;
    h1->Fill(tof, pin1e);
    if(fevent->LGBSize()>0){
      for(int i=0;i<fevent->LGFSize();i++){
        double lgft = fevent->LGF().at(i).GetTimestamp();
        double lgfe = fevent->LGF().at(i).GetEnergy();
        int    lgfn = fevent->LGF().at(i).GetNumber();
        for(int j=0;j<fevent->LGBSize();j++){
          double lgbt = fevent->LGB().at(j).GetTimestamp();
          double lgbe = fevent->LGB().at(j).GetEnergy();
          int    lgbn = fevent->LGB().at(j).GetNumber();
          if((lgfn>=64 && lgfn<=79) && lge[0]->IsInside(lgfe,lgbe)){
            dt = lgbt - lgft;
            t1->Fill(dt);
          }
          if((lgfn>=48 && lgfn<=63) && lge[1]->IsInside(lgfe,lgbe)){
            dt = lgbt - lgft;
            t1->Fill(dt);
          }
          if((lgfn>=40 && lgfn<=48) && lge[2]->IsInside(lgfe,lgbe)){
            dt = lgbt - lgft;
            t1->Fill(dt);
          }
        }
      }
    }
    if(dt!=0) h2->Fill(tof,pin1e);

    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
    
  }
  printf("on entry %lu / %lu \n",x,nentries);


  //TFile *outfile = new TFile(Form("histogram%i.root",runnum),"recreate");
  TFile *outfile = new TFile("beyondveto1_103.root","recreate");
  h1 ->Write();
  h2 ->Write();
  t1 ->Write();
  outfile->Close();


}




void quick(){

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate.root");
  TCutG *ban11 = (TCutG *)cutf2->Get("ban11");
  TCutG *ban12 = (TCutG *)cutf2->Get("ban12");
  TCutG *ban2  = (TCutG *)cutf2->Get("ban2");
  TCutG *ban31 = (TCutG *)cutf2->Get("ban31");
  TCutG *ban32 = (TCutG *)cutf2->Get("ban32");
 
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");


  int runnum[9] = {1030,1031,1032,1033,1034,1035,1036,1037,1039}; 
  TChain *chan = new TChain("event");
  for(int i=0;i<9;i++){
    chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i*.root",runnum[i]));
  }
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  TH1D *tdif1 = new TH1D("tdif1","t dif between 2 events", 1000,0,0);
  TH1D *tdif2 = new TH1D("tdif2","t dif between 2 events with lgf or lgb", 1000,0,0);
  TH1D *tdif3 = new TH1D("tdif3","t dif between 2 events with lgf", 1000,0,0);
  
  double time1 = -1;
  double time2 = -1;
  double time3 = -1;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    double dt = 0;
    // check veto
    bool veto = false;
    if(fevent->LGFSize()>0){
      double lgfmax = fevent->LGFMax().GetEnergy();
      if(fevent->SSSDSize()>0){
        for(int i=0;i<fevent->SSSDSize();i++){
          double sssdn = fevent->SSSD().at(i).GetStrip()-1;
          double sssdc = fevent->SSSD().at(i).GetCharge();
          double sssdt = fevent->SSSD().at(i).GetTimestamp();
          if(cutstrip->IsInside(sssdn, fevent->LGFMax().GetStrip()-1)){
            if(fevent->LGFMax1().GetEnergy()>0){
              if(ban11->IsInside(sssdc, lgfmax) || ban12->IsInside(sssdc, lgfmax)){
                veto = true;
                break;
              }
            } 
            if(fevent->LGFMax2().GetEnergy()>0){
              if(ban2->IsInside(sssdc, lgfmax)){
                veto = true;
                break;
              }
            } 
            if(fevent->LGFMax3().GetEnergy()>0){
              if(ban31->IsInside(sssdc, lgfmax) || ban32->IsInside(sssdc, lgfmax)){
                veto = true;
                break;
              }
            } 
          }
        }
      }
    }// end checking veto


    if(veto) continue;
    if(time1<0){
      time1 = fevent->Pin1T();
    }else{
      dt = fevent->Pin1T() - time1;
      dt = dt/1000000.;
      tdif1->Fill(dt);
      time1 = fevent->Pin1T(); 
    }
    if(fevent->LGFSize()>0 || fevent->LGBSize()>0){
      if(time2<0){
        time2 = fevent->Pin1T();
      }else{
        dt = fevent->Pin1T() - time2;
        dt = dt/1000000.;
        tdif2->Fill(dt);
        time2 = fevent->Pin1T(); 
      }
    }
    if(fevent->LGFSize()>0){
      if(time3<0){
        time3 = fevent->Pin1T();
      }else{
        dt = fevent->Pin1T() - time3;
        dt = dt/1000000.;
        tdif3->Fill(dt);
        time3 = fevent->Pin1T(); 
      }
    }



    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
    
  }
  printf("on entry %lu / %lu \n",x,nentries);


  //TFile *outfile = new TFile(Form("histogram%i.root",runnum),"recreate");
  TFile *outfile = new TFile("quick_103.root","recreate");
  tdif1 ->Write();
  tdif2 ->Write();
  tdif3 ->Write();
  outfile->Close();


}






