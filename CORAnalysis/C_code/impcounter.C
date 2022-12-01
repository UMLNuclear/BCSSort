

void counter(){


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

  long count0 = 0;
  long count1 = 0;
  long count2 = 0;
  long count3 = 0;
  long count4 = 0;
  long count5 = 0;
  long count6 = 0;
  long count7 = 0;


  double time = -1;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->LGFSize()==0){
      count1++;
      if(fevent->LGBSize()==0){
        count2++;
      }
    }
    if(fevent->LGBSize()==0) count6++;
    if(fevent->LGFSize()>2) count4++;
    if(fevent->LGBSize()>2) count3++;
    if(fevent->LGFSize()>0 && fevent->Pin2E()<0) count5++;


    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
  }
  printf("on entry %lu / %lu \n",x,nentries);
  
  std::cout << "total entries    = " << nentries << std::endl;
  std::cout << "implant          = " << nentries   << std::endl;
  std::cout << "(imp)no lgf      = " << count1 << std::endl;
  std::cout << "(imp)no lgb      = " << count6 << std::endl;
  std::cout << "(imp)no lg       = " << count2 << std::endl;
  std::cout << "(imp)lgfsize(>2) = " << count4 << std::endl;
  std::cout << "(imp)lgbsize(>2) = " << count3 << std::endl;
  std::cout << "(imp)no pin2     = " << count5 << std::endl;

}





void beyondvetocounter(){
  
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  TCutG *blob[6];
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }

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

  TChain *chan = new TChain("event");
  //int runnum[9] = {1040,1041,1042,1043,1044,1046,1047,1048,1049};
  //for(int i=0;i<9;i++){
  //  chan->Add(Form("/home/zhu/packages/BCSSort/junk/data_event/event%i*.root",runnum[i]));
  //}
  chan->Add("/home/zhu/packages/BCSSort/junk/data_event/event1040-00.root");
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  int count0[6];
  int count1[6];
  int count2[6];
  int count3[6];
  int count4[6];
  int count5[6];
  
  for(int i=0;i<6;i++){
    count0[i] = 0;
    count1[i] = 0;
    count2[i] = 0;
  }

  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->LGFSize()==0) continue;
    double pin1e = fevent->Pin1E();
    double tof   = fevent->I2S();
    double lgfmax = fevent->LGFMax().GetEnergy();
    double dt = 0;
    bool veto  = false;
    bool veto2 = false;
    for(int i=0;i<6;i++){
      if(blob[i]->IsInside(tof, pin1e)){
        count0[i]++;
      }
    }

    //check veto
    if(fevent->SSSDSize()>0){
      for(int i=0;i<fevent->SSSDSize();i++){
        double sssdn = fevent->SSSD().at(i).GetStrip()-1;
        double sssdc = fevent->SSSD().at(i).GetCharge();
        double sssdt = fevent->SSSD().at(i).GetTimestamp();
        if(fevent->LGFMax1().GetEnergy()>0){
          if(ban11->IsInside(sssdc, lgfmax) || ban12->IsInside(sssdc, lgfmax)){
            veto2 = true;
          }
        }
        if(fevent->LGFMax2().GetEnergy()>0){
          if(ban2->IsInside(sssdc, lgfmax)){
            veto2 = true;
          }
        }
        if(fevent->LGFMax3().GetEnergy()>0){
          if(ban31->IsInside(sssdc, lgfmax) || ban32->IsInside(sssdc, lgfmax)){
            veto2 = true;
          }
        }
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
    }//end checking veto
    
    if(!veto2){
      for(int i=0;i<6;i++){
        if(blob[i]->IsInside(tof, pin1e)){
          count1[i]++;
        }
      }
    }
    if(!veto){
      for(int i=0;i<6;i++){
        if(blob[i]->IsInside(tof, pin1e)){
          count2[i]++;
        }
      }
    }
    if((x%5000)==0){
      printf(" on entry = %lu / %lu \r", x, nentries);
      fflush(stdout);
    }
  

  }
  printf(" on entry = %lu / %lu \n", x, nentries);

  for(int i=0;i<6;i++){
    double dif = count2[i]-count1[i];
    dif = dif/count2[i]*100;
    printf("counts in %s: \n", blob[i]->GetName());
    printf("no gate = %i \n", count0[i]);
    printf("veto1   = %i \n", count2[i]);
    printf("veto2   = %i(miss %.2f%%) \n", count1[i], dif);
  } 



}


void quick(){

  TChain *chan = new TChain("event");
  //int runnum[9] = {1030,1031,1032,1033,1034,1035,1036,1037,1039};
  //for(int i=0;i<9;i++){
  //  chan->Add(Form("/home/zhu/packages/BCSSort/junk/banana_event/event%i*.root",runnum[i]));
  //}
  chan->Add("/home/zhu/packages/BCSSort/junk/event/example/veto/event1040-00.root");
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;
  
  int count0 = 0;

  TH1D *tdif = new TH1D("dt","dt",2000,0,20000);
  double time1 = -1;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    if(fevent->Pin1E()>0) count0++;
    if(x==0){
      time1 = fevent->fHits.at(fevent->Size()-1).GetTimestamp();
    }else{
      double dt = fevent->fHits.at(0).GetTimestamp() - time1;
      time1 = fevent->fHits.at(fevent->Size()-1).GetTimestamp();
      tdif->Fill(dt); 
    }
     
    if((x%5000)==0){
      printf("on entry = %lu / %lu \r",x,nentries);
      fflush(stdout);
    }  

  }

  printf("on entry %lu / %lu \n",x,nentries);

  cout<<"imp = "<<count0<<std::endl;

}


















