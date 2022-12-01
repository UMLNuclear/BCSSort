



void hist_banana2(){

  TChannel::ReadDetMapFile();
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate2.root");
  TCutG *ban = (TCutG *)cutf2->Get("ban");

  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  TCutG *cutstrip = (TCutG *)cutf3->Get("cutstrip");

  TH2D *h0 = 0;  
  TH2D *h1 = 0;  

  for(int runnum=1040;runnum<1050;runnum++){

    if(!h0){ h0 = new TH2D("pid",    "pid LGF>0 && LGB>0" , 1000,8000,18000,300,4200,7200); }
    if(!h1){ h1 = new TH2D("pid1",   "pid out of veto"    , 1000,8000,18000,300,4200,7200); }

    TChain *chan = new TChain("event");
    chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/implant/event%i*.root",runnum));
    long nentries = chan->GetEntries();
    long x = 0;
    BCSEvent *fevent = new BCSEvent;
    chan->SetBranchAddress("BCSEvent", &fevent);

    for(x=0;x<nentries;x++){
      fevent->Clear();
      chan->GetEntry(x);
      double pin1e = fevent->Pin1E();
      double tof   = fevent->I2S();     
      bool flag_veto = false;
      if(fevent->LGFSize()==0 || fevent->LGBSize()==0) continue;
      if(fevent->LGFSize()>10 || fevent->LGBSize()>10) continue;
      h0->Fill(tof, pin1e);
      
      int    lgfmaxn = fevent->LGFMax().GetStrip()-1;
      double lgfmaxe = fevent->LGFMax().GetEnergy();
      //veto check
      if(fevent->SSSDSize()>0){
        for(int i=0;i<fevent->SSSDSize();i++){
          int    sssdn = fevent->SSSD()[i].GetStrip()-1;
          double sssdc = fevent->SSSD()[i].GetCharge();
          if(cutstrip->IsInside(sssdn, lgfmaxn)){
            if(ban->IsInside(sssdc, lgfmaxe)){
              flag_veto = true;
              break;
            }
          }
        }
      }// veto check done

      //out of veto
      if(flag_veto) continue;
      h1->Fill(tof,pin1e);

      if((x%5000)==0){
        printf("run%i.root: on entry %lu / %lu \r",runnum,x,nentries);
        fflush(stdout);
      }

    }
    printf("run%i.root: on entry %lu / %lu \n",runnum,x,nentries);

  }







}
