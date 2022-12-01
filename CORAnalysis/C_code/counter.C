

void counter(){


  int runnum[3] = {1030,1031,1032};
  TChain *chan = new TChain("event");
  for(int i=0;i<3;i++){
    chan->Add(Form("/home/zhu/packages/BCSSort/junk/banana_event/event%i*.root",runnum[i]));
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

  TH2D *h1 = new TH2D("pid","pid", 1000,8000,18000,300,4200,7200);
  TH1D *h2 = new TH1D("hgfsize0", "hgfsize in dec",40,0,40);
  TH1D *h3 = new TH1D("hgbsize0", "hgbsize in dec",40,0,40);
  TH1D *h4 = new TH1D("hgfsize", "hgfsize with ge + hg",40,0,40);
  TH1D *h5 = new TH1D("hgbsize", "hgbsize with ge + hg",40,0,40);
  TH1D *h6 = new TH1D("dt1", "(dec)dt hgf - hgb",1000,-5000,5000);
  
  TH1D *g0 = new TH1D("size0","size of decay"                 , 1000,0,1000);
  TH1D *g1 = new TH1D("size1","size of decay(hg + ge + labr)" , 1000,0,1000);
  TH1D *g2 = new TH1D("size2","size of decay(hg + gam)"       , 1000,0,1000);
  TH1D *g3 = new TH1D("size3","size of decay(hg + no ge)"     , 1000,0,1000);
  TH1D *g4 = new TH1D("size4","size of decay(hg + no gam)"    , 1000,0,1000);
  TH1D *g5 = new TH1D("size5","size of decay(no hg + ge)"     , 1000,0,1000);
  TH1D *g6 = new TH1D("size6","size of decay(no hg + no gam)" , 1000,0,1000);
  TH1D *g7 = new TH1D("size7","size of decay(ge)"             , 1000,0,1000);

  double time = -1;
  int size0 = 0;
  int size1 = 0;
  int size2 = 0;
  int size3 = 0;
  int size4 = 0;
  int size5 = 0;
  int size6 = 0;
  int size7 = 0;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    //imp
    if(fevent->Pin1E()>0){
    }else{
      count0++;
      size0++;
      h2->Fill(fevent->HGFSize());
      h3->Fill(fevent->HGBSize());
      if((fevent->HGFSize()==0 && fevent->HGBSize()==0) && (fevent->HPGeSize()==0 && fevent->LaBrSize()==0)){
        count1++;
      }
      if(fevent->HPGeSize()>0){
        count2++;
        if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
          count3++;
          h4->Fill(fevent->HGFSize());
          h5->Fill(fevent->HGBSize());
        }
      }else{
        if(fevent->LaBrSize()>0){
          count4++;
        }
      }
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        count5++;
        if(fevent->HPGeSize()>0){
          count6++;
        }
        for(int i=0;i<fevent->HGFSize();i++){
          double hgft = fevent->HGF()[i].GetTimestamp();
          for(int j=0;j<fevent->HGBSize();j++){
            double hgbt = fevent->HGB()[j].GetTimestamp();
            double dt = hgft - hgbt;
            h6->Fill(dt);
          }
        }
      }     
    }






    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
  }
  printf("on entry %lu / %lu \n",x,nentries);
  
  std::cout << "total entries             = " << nentries << std::endl;
  std::cout << "decay                     = " << count0   << std::endl;
  std::cout << "(dec) no hg + no gam      = " << count1 << std::endl;
  std::cout << "(dec)ge                   = " << count2 << std::endl;
  std::cout << "(ge) hgf + hgb            = " << count3 << std::endl;
  std::cout << "(!ge)LaBr                 = " << count4 << std::endl;
  std::cout << "(dec)hgf + hgb            = " << count5 << std::endl;
  std::cout << "(hgf+hgb)ge               = " << count6 << std::endl;
}


void counter_corlgeo(){


  TChain *chan = new TChain("event");
  //int runnum[3] = {1030,1031,1032};
  //for(int i=0;i<3;i++){
  //  chan->Add(Form("/home/zhu/packages/BCSSort/junk/banana_event/event%i*.root",runnum[i]));
  //}
  chan->Add(Form("/home/zhu/packages/BCSSort/junk/banana_event/event%i-05.root",1032));
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

  TH1D *h1 = new TH1D("dt1", "dt1", 500,0,5000);

  pixel piximp = std::make_pair(-1,-1);
  
  //nentries = 6e4;
  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    //imp
    if(fevent->Pin1E()>0){
      count0++;
      piximp = fevent->LGPixel();
      if(piximp.first<0 || piximp.second<0){
        count1++;
        continue;
      }
      count5++;
      count6++;
    }else{
      if(piximp.first<0 || piximp.second<0) continue;
      count2++;
      bool flag = false;
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        count3++;
        for(int i=0;i<fevent->HGFSize();i++){
          int    hgfn = fevent->HGF()[i].GetStrip()-1;
          double hgft = fevent->HGF()[i].GetTimestamp();
          for(int j=0;j<fevent->HGBSize();j++){
            int    hgbn = fevent->HGB()[j].GetStrip()-1;
            double hgbt = fevent->HGB()[j].GetTimestamp();
            double dt = hgft- hgbt;
            if(dt>=30 && dt<140){
              h1->Fill(dt);
              if(fabs(hgfn-piximp.first)<=1 && fabs(hgbn-piximp.second)<=1){
                flag = true;
              }
            } 
          }
        }
      }
      if(flag){
        count4++;
        count6++;
      }

    }






    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,nentries);
      fflush(stdout);
    } 
  }
  printf("on entry %lu / %lu \n",x,nentries);
  
  std::cout << "total entries             = " << nentries << std::endl;
  std::cout << "implant                   = " << count0   << std::endl;
  std::cout << "(imp)no lgb               = " << count1   << std::endl;
  std::cout << "(imp)good                 = " << count5   << std::endl;
  std::cout << "decay                     = " << count2   << std::endl;
  std::cout << "(dec)hgf + hgb            = " << count3   << std::endl;
  std::cout << "(hgf+hgb)good correlation = " << count4   << std::endl;
  std::cout << "new tree entries          = " << count6   << std::endl;
}

void print(){
  int runnum = 1032;
  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/junk/banana_event/event%i-00*.root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  TChannel::ReadDetMapFile();
  long nentries = chan->GetEntries();
  long x = 0;

  nentries = 6e3;
  for(x=590;x<nentries;x++){
    fevent->Clear(); 
    chan->GetEntry(x);
    if(fevent->Pin1E()>0){
      printf("===============IMP===========\n");
      printf("entry = %lu\t LGFSize = %i \t LGBSize = %i\n", x, fevent->LGFSize(), fevent->LGBSize());
      if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
        for(int i=0;i<fevent->LGFSize();i++){
          int    lgfn = fevent->LGF()[i].GetStrip()-1;
          double lgfe = fevent->LGF()[i].GetEnergy();
          for(int j=0;j<fevent->LGBSize();j++){
            int    lgbn = fevent->LGB()[j].GetStrip()-1;
            double lgbe = fevent->LGB()[j].GetEnergy();
            printf("LG[%i,%i] = (%f, %f)\n", lgfn, lgbn, lgfe, lgbe);
          }
        }
      }
    }else{
      if(fevent->HGFSize()>0 && fevent->HGBSize()>0){
        //printf("entry = %lu\t HGFSize = %i \t HGBSize = %i\n", x, fevent->HGFSize(), fevent->HGBSize());
        bool flag1 = false;
        bool flag2 = false;
        for(int i=0;i<fevent->HGFSize();i++){
          int    hgfn = fevent->HGF()[i].GetStrip()-1;
          double hgfe = fevent->HGF()[i].GetEnergy();
          if(hgfn>=25 && hgfn<=27) {
            flag1 = true;
            break;
          }
        }
        for(int j=0;j<fevent->HGBSize();j++){
          int    hgbn = fevent->HGB()[j].GetStrip()-1;
          double hgbe = fevent->HGB()[j].GetEnergy();
          if(hgbn>=17 && hgbn<=21) {
            flag2 = true;
            break;
          } 
        }
        if(flag1 && flag2){
          for(int i=0;i<fevent->HGFSize();i++){
            int    hgfn = fevent->HGF()[i].GetStrip()-1;
            double hgfe = fevent->HGF()[i].GetEnergy();
            double hgft = fevent->HGF()[i].GetTimestamp() - fevent->Pin1T();
            for(int j=0;j<fevent->HGBSize();j++){
              int    hgbn = fevent->HGB()[j].GetStrip()-1;
              double hgbe = fevent->HGB()[j].GetEnergy();
              printf("entry = %lu(%f) \t HG[%i,%i] = (%.2f, %.2f) GeSize = %i\n", x, hgft, hgfn, hgbn, hgfe, hgbe, fevent->HPGeSize());
            }
          }
        }
      }
    }
  }



}







void quickcheck(){


  TChannel::ReadDetMapFile();
  for(int runnum=1045;runnum<1046;runnum++){
    for(int subrun=0;subrun<20;subrun++){
      TChain *chan = new TChain("event");
      chan->Add(Form("/home/zhu/packages/BCSSort/junk/event/example/veto/event%i-%02i.root",runnum, subrun));
      long nentries = chan->GetEntries();
      if(nentries==0){
        printf("event%i-%02i.root doesn't exist\n", runnum, subrun);
        continue;
      }
      long x = 0;
      BCSEvent *fevent = new BCSEvent;
      chan->SetBranchAddress("BCSEvent", &fevent);
      for(x=0;x<nentries;x++){
        fevent->Clear();
        chan->GetEntry(x);
        if(fevent->LGFSize()>3){
          printf("event%i-%02i: entry = %lu, LGFSize = %i\n", runnum, subrun, x, fevent->LGFSize());
        }
        if((x%5000)==0){
          printf("event%i-%02i: on entry %lu /%lu \r", runnum, subrun, x, nentries);
          fflush(stdout);
        }
  
      } 
      printf("event%i-%02i: on entry %lu /%lu \n", runnum, subrun, x, nentries);
    }
  }







}






