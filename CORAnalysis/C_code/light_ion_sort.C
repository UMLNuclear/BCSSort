
void hist(int runnum = 1040){
  TChannel::ReadDetMapFile();
  
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cuts_light.root");
  TCutG *lblob[6];
  for(int i=0;i<6;i++){
    lblob[i] = (TCutG *)cutf1->Get(Form("lblob%i",i));
  }

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf2->Get(Form("blob%i",i));
  }
  

  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/decay_hgcc_cut.root");
  std::map<std::pair<int,int>, TCutG *>hgcut;
  for(int i=1;i<4;i++){
    for(int j=1;j<4;j++){
      std::pair<int,int> cutpix = std::make_pair(i,j);
      hgcut[cutpix] = (TCutG*)cutf4->Get(Form("cut%i%i",i,j));
    }
  }
 

  TH2D *pid0 = new TH2D("pid","pid",2800,0,28000,1400,0,14000);
  TH2D *pidl[6];  // pid for light ions;
  TH2D *pidh[6];  // pid for heavy ions;
  TH2D *hmlgh[6]; // hitmap for lg for heavy ions
  TH2D *hmlgl[6]; // hitmap for lg for light ions
  TH2D *hmhgh[6]; // hitmap for hg heavy ions;
  TH2D *hmhgl0[6];// hitmap for hg light ions without lg fired;
  TH2D *hmhgl1[6];// hitmap for hg light ions with lg fired;
  TH2D *cc0l1[9]; //2D cc for light ion blob0 with lg fired
  TH2D *cc0l0[9]; //2D cc for light ion blob0 with lg no fired
  TH2D *cc1l1[9]; //2D cc for light ion blob1 with lg fired
  TH2D *cc1l0[9]; //2D cc for light ion blob1 with lg no fired
  TH2D *cc2l1[9]; //2D cc for light ion blob2 with lg fired
  TH2D *cc2l0[9]; //2D cc for light ion blob2 with lg no fired
  TH2D *cc3l1[9]; //2D cc for light ion blob3 with lg fired
  TH2D *cc3l0[9]; //2D cc for light ion blob3 with lg no fired
  TH2D *cc4l1[9]; //2D cc for light ion blob4 with lg fired
  TH2D *cc4l0[9]; //2D cc for light ion blob4 with lg no fired
  TH2D *cc5l1[9]; //2D cc for light ion blob5 with lg fired
  TH2D *cc5l0[9]; //2D cc for light ion blob5 with lg no fired

  for(int i=0;i<6;i++){
    pidl[i]  = new TH2D(Form("pidl%i",i), Form("pid light ion blob%i",i),2800,0,28000,1400,0,14000);  
    pidh[i]  = new TH2D(Form("pidh%i",i), Form("pid heavy ion blob%i",i),2800,0,28000,1400,0,14000);
    hmlgh[i] = new TH2D(Form("hmlg_h%i",i), Form("hitmap for lowgain, heavy ions blob%i",i), 40,0,40,40,0,40);  
    hmlgl[i] = new TH2D(Form("hmlg_l%i",i), Form("hitmap for lowgain, light ions blob%i",i), 40,0,40,40,0,40);  
    hmhgh[i] = new TH2D(Form("hmhg_h%i",i), Form("hitmap for highgain,with lowgain, heavy ions blob%i",i), 40,0,40,40,0,40);  
    hmhgl0[i]= new TH2D(Form("hmhg_nolg_l%i",i), Form("hitmap for highgain, no lowgain, light ions blob%i",i), 40,0,40,40,0,40);  
    hmhgl1[i]= new TH2D(Form("hmhg_lg_l%i",i), Form("hitmap for highgain, with lowgain, light ions blob%i",i), 40,0,40,40,0,40);  
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      cc0l1[3*i+j] = new TH2D(Form("cc0_lg%i%i",i+1,j+1),Form("light ions[0], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc0l0[3*i+j] = new TH2D(Form("cc0_nolg%i%i",i+1,j+1),Form("light ions[0], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc1l1[3*i+j] = new TH2D(Form("cc1_lg%i%i",i+1,j+1),Form("light ions[1], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc1l0[3*i+j] = new TH2D(Form("cc1_nolg%i%i",i+1,j+1),Form("light ions[1], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc2l1[3*i+j] = new TH2D(Form("cc2_lg%i%i",i+1,j+1),Form("light ions[2], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc2l0[3*i+j] = new TH2D(Form("cc2_nolg%i%i",i+1,j+1),Form("light ions[2], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc3l1[3*i+j] = new TH2D(Form("cc3_lg%i%i",i+1,j+1),Form("light ions[3], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc3l0[3*i+j] = new TH2D(Form("cc3_nolg%i%i",i+1,j+1),Form("light ions[3], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc4l1[3*i+j] = new TH2D(Form("cc4_lg%i%i",i+1,j+1),Form("light ions[4], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc4l0[3*i+j] = new TH2D(Form("cc4_nolg%i%i",i+1,j+1),Form("light ions[4], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc5l1[3*i+j] = new TH2D(Form("cc5_lg%i%i",i+1,j+1),Form("light ions[5], with lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
      cc5l0[3*i+j] = new TH2D(Form("cc5_nolg%i%i",i+1,j+1),Form("light ions[5], no lowgain, hgfc(%i) vs hgbc(%i)",i+1,j+1), 500,0,25000,500,0,25000);
    }
  } 

  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/allimp/event%i*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long x=0;
  long nentries = chan->GetEntries();
  
  //nentries = 1e5; 
  std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);
  for(x;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    bool flag_lg = false;
    if(fevent->HGFSize()==0 || fevent->HGBSize()==0) continue;
    if(fevent->LGFSize()>0 && fevent->LGBSize()>0) flag_lg = true;
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0];
    pid0->Fill(tof,pin1e);
    bool flag_l[6];
    bool flag_h[6];
    for(int m=0;m<6;m++){
      flag_l[m] = false;
      flag_h[m] = false;
    }
    for(int m=0;m<6;m++){
      if(blob[m] ->IsInside(tof,pin1e)) {
        flag_h[m] = true;  
        pidh[m]->Fill(tof,pin1e);
      }
      if(lblob[m]->IsInside(tof,pin1e)) {
        flag_l[m] = true; 
        pidl[m]->Fill(tof,pin1e);
      }
    }
    for(int i=0;i<fevent->HGFSize();i++){
      int    hgfn = fevent->HGF()[i].GetStrip()-1;
      double hgfc = fevent->HGF()[i].GetCharge();
      double hgft = fevent->HGF()[i].GetTimestamp();
      if(hgfc<10 || hgfc>25000) continue;
      for(int j=0;j<fevent->HGBSize();j++){
        int    hgbn = fevent->HGB()[j].GetStrip()-1;
        double hgbc = fevent->HGB()[j].GetCharge();
        double hgbt = fevent->HGB()[j].GetTimestamp();
        int car = 3*(hgfn/16)+hgbn/16;
        if(hgbc<10 || hgbc>25000) continue;
        std::pair<int,int> temp = std::make_pair((hgfn/16+1), (hgbn/16+1));
        double dt = hgft - hgbt;
        if(dt>=30 && dt<140){
          if(flag_l[0]) {
            if(flag_lg){hmhgl1[0]->Fill(hgfn,hgbn); cc0l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[0]->Fill(hgfn,hgbn); cc0l0[car]->Fill(hgfc,hgbc);}
          }
          if(flag_l[1]) {
            if(flag_lg){hmhgl1[1]->Fill(hgfn,hgbn); cc1l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[1]->Fill(hgfn,hgbn); cc1l0[car]->Fill(hgfc,hgbc);}
          }
          if(flag_l[2]) {
            if(flag_lg){hmhgl1[2]->Fill(hgfn,hgbn); cc2l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[2]->Fill(hgfn,hgbn); cc2l0[car]->Fill(hgfc,hgbc);}
          }
          if(flag_l[3]) {
            if(flag_lg){hmhgl1[3]->Fill(hgfn,hgbn); cc3l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[3]->Fill(hgfn,hgbn); cc3l0[car]->Fill(hgfc,hgbc);}
          }
          if(flag_l[4]) {
            if(flag_lg){hmhgl1[4]->Fill(hgfn,hgbn); cc4l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[4]->Fill(hgfn,hgbn); cc4l0[car]->Fill(hgfc,hgbc);}
          }
          if(flag_l[5]) {
            if(flag_lg){hmhgl1[5]->Fill(hgfn,hgbn); cc5l1[car]->Fill(hgfc,hgbc);}
            else       {hmhgl0[5]->Fill(hgfn,hgbn); cc5l0[car]->Fill(hgfc,hgbc);}
          }
          for(int m=0;m<6;m++){
            if(flag_h[i]) hmhgh[m]->Fill(hgfc,hgbc);
          }
        } // prompt t gate
      } // end HGB loop
    } // end HGF loop
    if(fevent->LGFMax().GetEnergy()<0 || fevent->LGBMax().GetEnergy()<0) continue;
    for(int i=0;i<fevent->LGFSize();i++){
      int lgfn = fevent->LGF()[i].GetStrip()-1;
      for(int j=0;j<fevent->LGBSize();j++){
        int lgbn = fevent->LGB()[j].GetStrip()-1;
        for(int m=0;m<6;m++){
          if(flag_l[m]) hmlgl[m]->Fill(lgfn,lgbn);
          if(flag_h[m]) hmlgh[m]->Fill(lgfn,lgbn);
        }
      } // end LGF loop
    } // end LGB loop

    if((x%5000)==0){
      printf("run%i.root; %lu / %lu\r",runnum,x,nentries);
      fflush(stdout);
    }
  }
  printf("run%i.root; %lu / %lu\n",runnum,x,nentries);
  TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/allimp/hist%i.root",runnum),"recreate");
  pid0->Write();
  for(int i=0;i<6;i++){
    pidl[i]  -> Write();
    pidh[i]  -> Write();
    hmlgh[i] -> Write();
    hmlgl[i] -> Write();
    hmhgh[i] -> Write();
    hmhgl0[i]-> Write();
    hmhgl1[i]-> Write();
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      cc0l1[3*i+j] -> Write();
      cc0l0[3*i+j] -> Write();
      cc1l1[3*i+j] -> Write();
      cc1l0[3*i+j] -> Write();
      cc2l1[3*i+j] -> Write();
      cc2l0[3*i+j] -> Write();
      cc3l1[3*i+j] -> Write();
      cc3l0[3*i+j] -> Write();
      cc4l1[3*i+j] -> Write();
      cc4l0[3*i+j] -> Write();
      cc5l1[3*i+j] -> Write();
      cc5l0[3*i+j] -> Write();
    }
  } 
  newf->Close();
  
}



void hist2(int runnum=1040){

  TChannel::ReadDetMapFile();
  
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cuts_light.root");
  TCutG *lblob[6];
  for(int i=0;i<6;i++){
    lblob[i] = (TCutG *)cutf1->Get(Form("lblob%i",i));
  }

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf2->Get(Form("blob%i",i));
  }
 
  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/hgfc_hbc_cuts.root");
  std::map<std::pair<int,int>, TCutG *> hgc0;
  std::map<std::pair<int,int>, TCutG *> hgc1;
  std::map<std::pair<int,int>, TCutG *> hgc2;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      std::pair<int,int> temp = std::make_pair(i+1,j+1);
      hgc0[temp] = (TCutG *)cutf3->Get(Form("cc%i%i_%i",i+1,j+1,0));
      hgc1[temp] = (TCutG *)cutf3->Get(Form("cc%i%i_%i",i+1,j+1,1));
      hgc2[temp] = (TCutG *)cutf3->Get(Form("cc%i%i_%i",i+1,j+1,2));
    }
  }
 

  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/decay_hgcc_cut.root");
  std::map<std::pair<int,int>, TCutG *>hgcut;
  for(int i=1;i<4;i++){
    for(int j=1;j<4;j++){
      std::pair<int,int> cutpix = std::make_pair(i,j);
      hgcut[cutpix] = (TCutG*)cutf4->Get(Form("cut%i%i",i,j));
    }
  }

  TH2D *cc0[9];
  TH2D *cc1[9];
  TH2D *cc2[9];
  TH2D *cc3[9];
  TH2D *cc4[9];
  TH2D *cc5[9];
  TH2D *ct0[9];
  TH2D *ct1[9];
  TH2D *ct2[9];
  TH2D *ct3[9];
  TH2D *ct4[9];
  TH2D *ct5[9];
  TH2D *hml[6];
  TH2D *hmt[6];
  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      cc0[3*i+j] = new TH2D(Form("cc0_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",0,i+1,j+1),500,0,25000,500,0,25000);
      cc1[3*i+j] = new TH2D(Form("cc1_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",1,i+1,j+1),500,0,25000,500,0,25000);
      cc2[3*i+j] = new TH2D(Form("cc2_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",2,i+1,j+1),500,0,25000,500,0,25000);
      cc3[3*i+j] = new TH2D(Form("cc3_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",3,i+1,j+1),500,0,25000,500,0,25000);
      cc4[3*i+j] = new TH2D(Form("cc4_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",4,i+1,j+1),500,0,25000,500,0,25000);
      cc5[3*i+j] = new TH2D(Form("cc5_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i)",5,i+1,j+1),500,0,25000,500,0,25000);
      ct0[3*i+j] = new TH2D(Form("ct0_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",0,i+1,j+1),500,0,25000,500,0,25000);
      ct1[3*i+j] = new TH2D(Form("ct1_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",1,i+1,j+1),500,0,25000,500,0,25000);
      ct2[3*i+j] = new TH2D(Form("ct2_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",2,i+1,j+1),500,0,25000,500,0,25000);
      ct3[3*i+j] = new TH2D(Form("ct3_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",3,i+1,j+1),500,0,25000,500,0,25000);
      ct4[3*i+j] = new TH2D(Form("ct4_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",4,i+1,j+1),500,0,25000,500,0,25000);
      ct5[3*i+j] = new TH2D(Form("ct5_%i%i",i+1,j+1),Form("lblob[%i]: hgfc(%i) vs hgbc(%i) time gated",5,i+1,j+1),500,0,25000,500,0,25000);
    }
  }

  for(int i=0;i<6;i++){
    hml[i] = new TH2D(Form("hml%i",i),Form("lblob[%i]: hitmap",i), 40,0,40,40,0,40);
    hmt[i] = new TH2D(Form("hmt%i",i),Form("lblob[%i]: hitmap time gated",i), 40,0,40,40,0,40);
  }


  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/allimp/event%i*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long x=0;
  long nentries = chan->GetEntries();

  std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);

  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->HGFSize()==0 || fevent->HGBSize()==0) continue;
    if(fevent->HGFSize(10,25000)>3  || fevent->HGBSize(10,25000)>3)  continue;
    if(fevent->LGFSize()>0  || fevent->LGBSize()>0)  continue;
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    tof = tofpar[2]*tof*tof + tofpar[1]*tof + tofpar[0];
    for(int i=0;i<fevent->HGFSize();i++){
      int    hgfn = fevent->HGF()[i].GetStrip()-1;
      double hgfc = fevent->HGF()[i].GetCharge();
      double hgft = fevent->HGF()[i].GetTimestamp();
      if(hgfc<10 || hgfc>25000) continue;
      for(int j=0;j<fevent->HGBSize();j++){
        int    hgbn = fevent->HGB()[j].GetStrip()-1;
        double hgbc = fevent->HGB()[j].GetCharge();
        double hgbt = fevent->HGB()[j].GetTimestamp();
        int car = 3*(hgfn/16)+hgbn/16;
        if(hgbc<10 || hgbc>25000) continue;
        std::pair<int,int> temp = std::make_pair((hgfn/16+1), (hgbn/16+1));
        double dt = hgft - hgbt;
        if(lblob[0]->IsInside(tof,pin1e)) { cc0[car]->Fill(hgfc,hgbc); hml[0]->Fill(hgfn,hgbn); }
        if(lblob[1]->IsInside(tof,pin1e)) { cc1[car]->Fill(hgfc,hgbc); hml[1]->Fill(hgfn,hgbn); }
        if(lblob[2]->IsInside(tof,pin1e)) { cc2[car]->Fill(hgfc,hgbc); hml[2]->Fill(hgfn,hgbn); }
        if(lblob[3]->IsInside(tof,pin1e)) { cc3[car]->Fill(hgfc,hgbc); hml[3]->Fill(hgfn,hgbn); }
        if(lblob[4]->IsInside(tof,pin1e)) { cc4[car]->Fill(hgfc,hgbc); hml[4]->Fill(hgfn,hgbn); }
        if(lblob[5]->IsInside(tof,pin1e)) { cc5[car]->Fill(hgfc,hgbc); hml[5]->Fill(hgfn,hgbn); }
        if(dt>=30 && dt<140){
          if(lblob[0]->IsInside(tof,pin1e)) { ct0[car]->Fill(hgfc,hgbc); hmt[0]->Fill(hgfn,hgbn);}
          if(lblob[1]->IsInside(tof,pin1e)) { ct1[car]->Fill(hgfc,hgbc); hmt[1]->Fill(hgfn,hgbn);}
          if(lblob[2]->IsInside(tof,pin1e)) { ct2[car]->Fill(hgfc,hgbc); hmt[2]->Fill(hgfn,hgbn);}
          if(lblob[3]->IsInside(tof,pin1e)) { ct3[car]->Fill(hgfc,hgbc); hmt[3]->Fill(hgfn,hgbn);}
          if(lblob[4]->IsInside(tof,pin1e)) { ct4[car]->Fill(hgfc,hgbc); hmt[4]->Fill(hgfn,hgbn);}
          if(lblob[5]->IsInside(tof,pin1e)) { ct5[car]->Fill(hgfc,hgbc); hmt[5]->Fill(hgfn,hgbn);}
        }
      }   
    }
    
    if((x%5000)==0){
      printf("run%i.root: on entry %lu / %lu \r", runnum, x, nentries);
      fflush(stdout);
    }
  
  }
  printf("run%i.root: on entry %lu / %lu \n", runnum, x, nentries);

  TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/allimp/StripSizeGated/hist%i.root",runnum),"recreate");
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      cc0[3*i+j] -> Write();
      cc1[3*i+j] -> Write();
      cc2[3*i+j] -> Write();
      cc3[3*i+j] -> Write();
      cc4[3*i+j] -> Write();
      cc5[3*i+j] -> Write();
      ct0[3*i+j] -> Write();
      ct1[3*i+j] -> Write();
      ct2[3*i+j] -> Write();
      ct3[3*i+j] -> Write();
      ct4[3*i+j] -> Write();
      ct5[3*i+j] -> Write();
    }
  }

  for(int i=0;i<6;i++){
    hml[i] -> Write();
    hmt[i] -> Write();
  }
  newf->Close();


}
