
void print(){
  TChannel::ReadDetMapFile();
  
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cuts_light.root");
  TCutG *lblob[6];
  for(int i=0;i<6;i++){
    lblob[i] = (TCutG *)cutf1->Get(Form("lblob%i",i));
  }

  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/decay_hgcc_cut.root");
  std::map<std::pair<int,int>, TCutG *>hgcut;
  for(int i=1;i<4;i++){
    for(int j=1;j<4;j++){
      std::pair<int,int> cutpix = std::make_pair(i,j);
      hgcut[cutpix] = (TCutG*)cutf4->Get(Form("cut%i%i",i,j));
    }
  }

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/hgfc_hbc_cuts.root");
  std::map<std::pair<int,int>, std::vector<TCutG *>> hgc;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      std::pair<int,int> temp = std::make_pair(i+1,j+1);
      for(int m=0;m<3;m++){
        hgc[temp].push_back((TCutG *)cutf2->Get(Form("cc%i%i_%i",i+1,j+1,m)));  
      }
    }
  }

  TChain *chan = new TChain("event");
  chan->Add("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/decay/highgain/event1040*.root");
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long x=0.998e6;
  //long nentries = chan->GetEntries();
  long nentries = 1e6;

  for(x;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    bool flag = false;
    for(int i=0;i<fevent->HGFSize();i++){
      int    hgfn = fevent->HGF()[i].GetStrip()-1;
      double hgfc = fevent->HGF()[i].GetCharge();
      double hgft = fevent->HGF()[i].GetTimestamp();
      if(hgfc<10 || hgfc>25000) continue;
      for(int j=0;j<fevent->HGBSize();j++){
        int    hgbn = fevent->HGB()[j].GetStrip()-1;
        double hgbc = fevent->HGB()[j].GetCharge();
        double hgbt = fevent->HGB()[j].GetTimestamp();
        if(hgbc<10 || hgbc>25000) continue;
        std::pair<int,int> temp = std::make_pair((hgfn/16+1), (hgbn/16+1));
        double dt = hgft - hgbt;
        if(dt>=30 && dt<140){
          for(int m=0;m<3;m++){
            if(hgc[temp][m]->IsInside(hgfc,hgbc)){
              if(!flag)printf("\nEntry \t c.zone  HGF.N \t HGB.N \t HGF.C \t HGB.C\n");
              printf("%lu \t %i \t %i \t %i \t %.1f \t %.1f\n", x, m, hgfn, hgbn, hgfc, hgbc);
              flag = true;
            }
          } // hgfc_hgbc gate
        } // prompt t gate
      } // end HGB loop
    } // end HGF loop
  }
}

void draw(TChain *chan, BCSEvent *fevent, long x = 0){
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cuts_light.root");
  TCutG *lblob[6];
  for(int i=0;i<6;i++){
    lblob[i] = (TCutG *)cutf1->Get(Form("lblob%i",i));
  }

  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/decay_hgcc_cut.root");
  std::map<std::pair<int,int>, TCutG *>hgcut;
  for(int i=1;i<4;i++){
    for(int j=1;j<4;j++){
      std::pair<int,int> cutpix = std::make_pair(i,j);
      hgcut[cutpix] = (TCutG*)cutf4->Get(Form("cut%i%i",i,j));
    }
  }

  TH2D *hg = 0;
  
  fevent->Clear();
  chan->GetEntry(x);
  for(int i=0;i<fevent->HGFSize();i++){
    int    hgfn = fevent->HGF()[i].GetStrip()-1;
    double hgfc = fevent->HGF()[i].GetCharge();
    double hgft = fevent->HGF()[i].GetTimestamp();
    if(hgfc<10 || hgfc>25000) continue;
    for(int j=0;j<fevent->HGBSize();j++){
      int    hgbn = fevent->HGB()[j].GetStrip()-1;
      double hgbc = fevent->HGB()[j].GetCharge();
      double hgbt = fevent->HGB()[j].GetTimestamp();
      if(hgbc<10 || hgbc>25000) continue;
      std::pair<int,int> temp = std::make_pair((hgfn/16+1), (hgbn/16+1));
      double dt = hgft - hgbt;
      if(dt>=30 && dt<140){
        if(hgcut[temp]->IsInside(hgfc,hgbc)){
          if(!hg) hg = new TH2D(Form("HG_%lu",x), Form("HG_%lu",x),80,0,80,40,0,40);
          hg->Fill(hgbn,hgfn,hgfc);
          hg->Fill(hgbn+40,hgfn,hgbc);
        } // hgfc_hgbc gate
      } // prompt t gate
    } // end HGB loop
  } // end HGF loop
  if(!hg) return;
  TCanvas *c = new TCanvas;
  c->Divide(1,2);
  c->cd(1);
  fevent->DrawLG()->Draw();
  c->cd(2);
  hg->Draw();
  c->Modified();
  c->Update();
}

