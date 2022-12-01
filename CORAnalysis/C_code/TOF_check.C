



/*void TOF_check(){

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut.root"); 
  TCutG *na = (TCutG *)cutf1->Get("Na");
  TCutG *ne = (TCutG *)cutf1->Get("Ne");

  TChannel::ReadDetMapFile();

  TH2D *h11 = 0;
  TH2D *h12 = 0;
  TH2D *h21 = 0;
  TH2D *h22 = 0;

  int runnum;
  int subrun;
  for(runnum = 1004;runnum<1040;runnum++){
    std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);
    if(tofpar.empty()){
      printf("run%i_%i doesn't have tof correction par\n",runnum,subrun);
      continue;
    }
    for(subrun = 0;subrun<20;subrun++){

      TChain *chan = new TChain("event");
      chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i-%02d.root",runnum,subrun));
      BCSEvent *fevent = new BCSEvent;
      chan->SetBranchAddress("BCSEvent", &fevent);
      long nentries = chan->GetEntries();
      long x = 0;
      if(nentries == 0){
        printf("run%i_%02i doesn't exist\n",runnum,subrun);
        continue;
      }

      h11 = new TH2D(Form("toft0na%i_%i",runnum,subrun),"Na running time vs tof before correction",4000,0,4e3,1000,8000,18000); 
      h12 = new TH2D(Form("toft1na%i_%i",runnum,subrun),"Na running time vs tof after correction", 4000,0,4e3,1000,8000,18000); 
      h21 = new TH2D(Form("toft0ne%i_%i",runnum,subrun),"Ne running time vs tof before correction",4000,0,4e3,1000,8000,18000); 
      h22 = new TH2D(Form("toft1ne%i_%i",runnum,subrun),"Ne running time vs tof after correction", 4000,0,4e3,1000,8000,18000); 

      int count0 = 0;
      int time =-1;

      for(x=0;x<nentries;x++){
        fevent->Clear();
        chan->GetEntry(x);

        if(fevent->I2S()<0){
          count0++;
        }else{
          double tof = fevent->I2S();
          tof = tofpar[0]+tofpar[1]*tof+tofpar[2]*tof*tof;
          double runtime = fevent->Pin1T() - time;
          runtime = runtime/1e9;
          if(na->IsInside(fevent->I2S(), fevent->Pin1E())){
            h11->Fill(runtime, fevent->I2S());
            h12->Fill(runtime, tof);
          }
          if(ne->IsInside(fevent->I2S(), fevent->Pin1E())){
            h21->Fill(runtime, fevent->I2S());
            h22->Fill(runtime, tof);
          }
        }

        if((x%5000)==0){
          printf("on entry %lu / %lu \r",x,nentries);
          fflush(stdout);
        }

      }
      printf("on entry %lu / %lu \n",x,nentries);
      if(count0>0){
        printf("run%i-%i has %i entries without I2S\n\n", runnum, subrun, count0);
      }


      TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/junk/histogram/event_tof_correction/tofcheck%i_%i.root",runnum,subrun),"recreate");
      h11->Write();
      h12->Write();
      h21->Write();
      h22->Write();
      newf->Close();
    }
  }

}*/

void TOF_check(){

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut.root"); 
  TCutG *na = (TCutG *)cutf1->Get("Na");
  TCutG *ne = (TCutG *)cutf1->Get("Ne");

  TChannel::ReadDetMapFile();
  std::map<int,double[4][4]> xtalmat = Histogram::Get()->ReadMat();  

  TH2D *h11 = 0;
  TH2D *h12 = 0;
  TH2D *h21 = 0;
  TH2D *h22 = 0;

  int runnum;
  for(runnum = 1023;runnum<1139;runnum++){
    std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);
    if(tofpar.empty()){
      printf("run%i doesn't have tof correction par\n",runnum);
      continue;
    }
    TChain *chan = new TChain("beta");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/correlation5/beta%i*.root",runnum));
    Beta *fbeta = new Beta;
    chan->SetBranchAddress("Beta", &fbeta);
    long nentries = chan->GetEntries();
    long x = 0;
    if(nentries == 0){
      printf("run%i doesn't exist\n",runnum);
      continue;
    }

    h21 = new TH2D(Form("toft0ne%i",runnum),"Ne running time vs tof gated 150keV",4000,0,4e3,1000,8000,18000); 
    h22 = new TH2D(Form("toft1ne%i",runnum),"Ne running time vs tof gated 150keVBG", 4000,0,4e3,1000,8000,18000); 

    h11 = new TH2D(Form("tofne%i",runnum),"Ne running time vs tof",4000,0,4e3,1000,8000,18000);
    int count0 = 0;
    int time =-1;

    for(x=0;x<nentries;x++){
      fbeta->Clear();
      chan->GetEntry(x);
      double tof  = fbeta->fImplant.fI2S;
      double pin1e= fbeta->fImplant.fPIN1E;
      double runtime = fbeta->fImplant.fPIN1T - time;
      runtime = runtime/1e9;
      if(ne->IsInside(tof,pin1e)) {
        h11->Fill(runtime,tof) ;
      }

      if((x%5000)==0){
        printf("runnum%i: on entry %lu / %lu \r",runnum,x,nentries);
        fflush(stdout);
      }

    }
    printf("runnum%i: on entry %lu / %lu \n",runnum,x,nentries);

    TFile *newf = new TFile(Form("/home/zhu/packages/BCSSort/CORAnalysis/histogram/betasort5/TOFCheck/tofcheck%i.root",runnum),"recreate");
    h11->Write();
    newf->Close();
  }

}

void TOF_check1(){

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut.root"); 
  TCutG *na = (TCutG *)cutf1->Get("Na");
  TCutG *ne = (TCutG *)cutf1->Get("Ne");

  TChannel::ReadDetMapFile();
  std::map<int,double[4][4]> xtalmat = Histogram::Get()->ReadMat();  

  TH2D *h1 = new TH2D("tofne1","Ne within 10ms runnum vs tof gated 373keV"   , 110,0,110,1000,8000,18000);
  TH2D *h2 = new TH2D("tofne2","Ne within 10ms runnum vs tof gated 373keVBGR", 110,0,110,1000,8000,18000);
  TH2D *h3 = new TH2D("tofne3","Ne within 30ms runnum vs tof gated 373keV"   , 110,0,110,1000,8000,18000);
  TH2D *h4 = new TH2D("tofne4","Ne within 30ms runnum vs tof gated 373keVBGR", 110,0,110,1000,8000,18000);

  int runnum;
  int count = -1;
  for(runnum = 1023;runnum<1139;runnum++){
    std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);
    if(tofpar.empty()){
      printf("run%i doesn't have tof correction par\n",runnum);
      continue;
    }
    TChain *chan = new TChain("beta");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/correlation5/beta%i*.root",runnum));
    Beta *fbeta = new Beta;
    chan->SetBranchAddress("Beta", &fbeta);
    long nentries = chan->GetEntries();
    long x = 0;
    if(nentries == 0){
      printf("run%i doesn't exist\n",runnum);
      continue;
    }

    count++;
    for(x=0;x<nentries;x++){
      fbeta->Clear();
      chan->GetEntry(x);
      double tof  = fbeta->fImplant.fI2S;
      double pin1e= fbeta->fImplant.fPIN1E;
      
      for(int d=0;d<fbeta->DecaySize();d++){
        double decayt = fbeta->fDecay[d].fDecayTime;
        if(decayt>30) continue;
        std::map<int, Clover> clmap;
        for(int m=0;m<fbeta->fDecay[d].GeSize();m++){
          double em = fbeta->fDecay[d].fGe[m].GetEnergy();
          if(em<10 || em>4000) continue;
          int clnum = (fbeta->fDecay[d].fGe[m].GetNumber()-208)/4;
          clmap[clnum].Add(fbeta->fDecay[d].fGe[m]);
        }//gamma loop end;
      
        std::map<int, Clover>::iterator it1;
        for(it1=clmap.begin();it1!=clmap.end();it1++){
          double adde1 = 0;
          for(int cc1=0;cc1<it1->second.Size();cc1++){
            adde1 += it1->second.fXtal[cc1].GetEnergy();
            for(int cc2=0;cc2<it1->second.Size();cc2++){
              int num1 = (it1->second.fXtal[cc1].GetNumber()-208)%4;
              int num2 = (it1->second.fXtal[cc2].GetNumber()-208)%4;
              adde1 += (it1->second.fXtal[cc2].GetEnergy())*xtalmat[it1->first][num1][num2];
            }
          }
          it1->second.SetAddE(adde1);
          double gamEm = it1->second.AddbackE();
          if(gamEm>370 && gamEm<376){
            if(ne->IsInside(tof,pin1e)){
              if(decayt<10) {h1->Fill(count, tof);}
              if(decayt<30) {h3->Fill(count, tof);}
            }
          }
          if(gamEm>376 && gamEm<382){
            if(ne->IsInside(tof,pin1e)){
              if(decayt<10) {h2->Fill(count, tof);}
              if(decayt<30) {h4->Fill(count, tof);}
            }
          }
        } //addback loop end;
      }//decay loop end

      if((x%5000)==0){
        printf("runnum%i: on entry %lu / %lu \r",runnum,x,nentries);
        fflush(stdout);
      }

    }
    printf("runnum%i: on entry %lu / %lu \n",runnum,x,nentries);

  }
    TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/histogram/betasort5/TOFCheck/tofcheck_gate373.root","recreate");
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    newf->Close();

}

void TOF_gammaGate(){

  TFile *cutf = new TFile("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *na = (TCutG *)cutf->Get("Na");
  TCutG *ne = (TCutG *)cutf->Get("Ne");
  na->SetName("na");
  ne->SetName("ne");

  TH2D *pid[110];
  TH2D *hbg[110];

  int i=0;
  int runnum = 1022;
  while(i<110){
    runnum ++;
    if(runnum==1025 || runnum==1030 || runnum==1051 || runnum==1052 || runnum==1104 || runnum==1115) continue;
    TFile *infile = TFile::Open(Form("/home/zhu/packages/BCSSort/CORAnalysis/histogram/betasort5/pid/betahist%i.root",runnum));
    pid[i] = (TH2D *)infile->Get("pid4_1"); //PID gated 373keV within 10ms;
    hbg[i] = (TH2D *)infile->Get("pid5_1"); //PID gated 373keV within 10ms;
    pid[i] -> SetName(Form("pid_%i",runnum));
    hbg[i] -> SetName(Form("hbg_%i",runnum));
    i++;
  }

  TH1D *h1[110];
  TH1D *bg[110];
  TH1D *c1[110];
  TH1D *c2[11];

  for(int i=0;i<110;i++){
    h1[i] = pid[i]->ProjectionX(Form("tof%i",i), 1,300, "[ne]");
    bg[i] = hbg[i]->ProjectionX(Form("bg%i",i), 1,300, "[ne]");
    h1[i] -> Rebin(8);
    bg[i] -> Rebin(8);
    c1[i] = (TH1D *)h1[i]->Clone(Form("c%s",h1[i]->GetName()));
    c1[i] -> Add(bg[i],-1);
    int j=i/10;
    if((i%10)==0){
      c2[j] = (TH1D *)c1[i]->Clone(Form("pid_c%i",j));
    }else{
      c2[j]->Add(c1[i],1);
    }
  }

  TCanvas *c = new TCanvas;
  c->Divide(4,3);
  for(int i=0;i<11;i++){
    c2[i]->GetXaxis()->SetRangeUser(10000,14000);
    c->cd(i+1);
    c2[i]->Draw("hist");
    c->Update();
  }

}

void TOF_Draw(){

  TH2D *pid[110];

  int i=0;
  int runnum = 1022;
  while(i<110){
    runnum ++;
    if(runnum==1025 || runnum==1030 || runnum==1051 || runnum==1052 || runnum==1104 || runnum==1115) continue;
    TFile *infile = TFile::Open(Form("/home/zhu/packages/BCSSort/CORAnalysis/histogram/betasort5/TOFCheck/tofcheck%i.root",runnum));
    pid[i] = (TH2D *)infile->Get(Form("tofne%i",runnum)); 
    i++;
  }
  
  TCanvas *c1 = new TCanvas;
  TCanvas *c2 = new TCanvas;
  TCanvas *c3 = new TCanvas;
  TCanvas *c4 = new TCanvas;
  TCanvas *c5 = new TCanvas;
  c1->Divide(6,4);  
  c2->Divide(6,4);  
  c3->Divide(6,4);  
  c4->Divide(6,4);  
  c5->Divide(6,4);  

  for(int i=0;i<110;i++){
    pid[i]->RebinY(8);
    int m = i/22;
    int n = i%22;
    if(m==0){
      c1->cd(n+1);
      pid[i]->Draw("colz");
    }
    if(m==1){
      c2->cd(n+1);
      pid[i]->Draw("colz");
    }
    if(m==2){
      c3->cd(n+1);
      pid[i]->Draw("colz");
    }
    if(m==3){
      c4->cd(n+1);
      pid[i]->Draw("colz");
    }
    if(m==4){
      c5->cd(n+1);
      pid[i]->Draw("colz");
    }
    c1->Update();
    c2->Update();
    c3->Update();
    c4->Update();
    c5->Update();
  }

  TH1D *h[110];
  for(int i=0;i<110;i++){
    h[i] = pid[i]->ProjectionY(Form("tof%i",i), 1,4000);
  }
  TCanvas *c = new TCanvas;
  c->Divide(4,3);
  for(int i=0;i<110;i++){
    int m = i/10;
    int n = i%10;
    h[i]->SetLineColor(n);
    //h[i]->Rebin(8);
    if(n==0){
      c->cd(m+1);
      h[i]->Draw();
    }else{
      h[i]->Draw("same");
    }
    c->Update();
  }
  


}





