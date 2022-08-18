



void loop(){
  TChannel::ReadDetMapFile(); 
  TChain *cheve;
  BCSEvent *fevent;

  TChain *chimp;
  Implant *fimp;

  std::vector<TCutG*> veccut;
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/prompt_Imp104.root");
  TCutG *tof_dt = (TCutG *)cutf->Get("tof_dt");
  TCutG *pin1E_dt = (TCutG *)cutf->Get("pin1E_dt");
  veccut.push_back(tof_dt); // TOF vs timedif for 44S;
  veccut.push_back(pin1E_dt); // Pin1E vs timedif for S;

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pidcut104.root");
  TCutG *Na = (TCutG *)cutf2->Get("Na");
  TCutG *Ne = (TCutG *)cutf2->Get("Ne");

  ofstream ofile;
  ofile.open("good_imp_counter.txt");
  ofile<<setw(5)<<"runnum"     <<"\t"
       <<setw(5)<<"2gates "    <<"\t"
       <<setw(5)<<"2gates+fbl "<<"\t"
       <<setw(5)<<"beta_imp "<<"\t"
       <<setw(5)<<"mul_pin1 "<<"\t"
       <<setw(5)<<"mul_tof "<<std::endl;

  for(int runnum=1023;runnum<1140;runnum++){
    cheve = new TChain("event");
    fevent = new BCSEvent;
    cheve->Add(Form("data/5us_tofcor/event/event%i*",runnum));
    cheve->SetBranchAddress("BCSEvent", &fevent);

    chimp = new TChain("implant");
    fimp = new Implant;
    chimp->Add(Form("data/5us_tofcor/implant/imp_good/implant%i*",runnum));
    chimp->SetBranchAddress("Implant", &fimp);

    long x=0;
    long neve=cheve->GetEntries();  
    long nimp=chimp->GetEntries();  

    if(neve<1 || nimp<1) continue;


    double low, high, dt=0;
    long c1 = 0;
    long c2 = 0;
    long c3 = 0;
    long c4 = 0;
    long c5 = 0;
    long c6 = 0;
    long c7 = 0;
    long c8 = 0;
    long c9 = 0;
    long c10 = 0;
    long c11 = 0;
    long c12 = 0;
    long c13 = 0;
    long c14 = 0;
    long c15 = 0;
    long c16 = 0;
    long c17 = 0;
    long c18 = 0;
    long c19 = 0;
    long c20 = 0;
    long c21 = 0;
    long c22 = 0;
    long c23 = 0;
    long c24 = 0;
    long c25 = 0;
    long c26 = 0;
    long c27 = 0;

    bool cflag = false;
    bool hasDSSD = false;


    for(x=0;x<neve;x++){
      cheve->GetEntry(x);
      if(fevent->Pin1E()<=0) continue;
      if(fevent->DSSDloT()<=0) continue;
      dt = fevent->DSSDloT() - fevent->Pin1T();
      dt = dt/1000.;
      if(veccut[0]->IsInside(dt,fevent->I2S()) && veccut[1]->IsInside(dt,fevent->Pin1E())){
        c1++;
        if(fevent->LGFSize()>0 && fevent->LGBSize()>0){
          c2++;
          c10 = 0;
          c11 = 0;
          for(auto &it:fevent->fHits){
            if(it.GetNumber()==181) c10++;
            if(it.GetNumber()==177) c11++;
          }
          if(c10>1) c5++;
          if(c11>1) c6++;
        }
      } 
    }

    for(x=0;x<nimp;x++){
      chimp->GetEntry(x);
      if(fimp->fDSSDFront.size()>0 && fimp->fDSSDBack.size()>0){
        c3++;
        dt = fimp->DSSDloT() - fimp->fPIN1T;
        dt = dt/1000.;
        if(veccut[0]->IsInside(dt,fimp->fI2S) && veccut[1]->IsInside(dt,fimp->fPIN1E)){
          c4++;
        }
      }  
    }



    //std::cout<<std::endl;
    //std::cout<< "(Eve) good imp gate*2= " << c1 << std::endl;
    //std::cout<< "(Eve) good imp gate*2 + LGF&&LGB>0 = " << c2 << std::endl;
    //std::cout<< "(Eve) good imp gate*2 + LGF&&LGB>0 + mul_pin1= " << c5 << std::endl;
    //std::cout<<std::endl;
    //std::cout<< "(Good Imp) total entries = " << nimp << std::endl;
    //std::cout<< "(Good Imp) LGF && LGB >0 = " << c3 << std::endl;
    //std::cout<< "(Good Imp) LGF && LGB >0 && FilterPixel = " << c4 << std::endl;
    //std::cout<<std::endl;

    ofile<<setw(5)<<runnum<<"\t\t"
         <<setw(5)<<c1    <<"\t\t" 
         <<setw(5)<<c3    <<"\t\t" 
         <<setw(5)<<c4    <<"\t\t"
         <<setw(5)<<c5    <<"\t\t"
         <<setw(5)<<c6    <<std::endl;

  }
  ofile.close();
  return;

}







