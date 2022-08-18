

#include<util.h>

void betaprint(){

  TChannel::ReadDetMapFile();

  ofstream ofile;
  ofile.open("beta2_GetPixel_check.txt");
  ofile << "run"       << "\t\t"
    << "imp_Ne"    << "\t\t"
    << "decay"     << "\t\t"
    << "imp_bad"   << "\t\t"
    << "decay_bad" << std::endl;
  
  TFile *mycut1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)mycut1->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");
  TChain *beta; 

  for(int runnum=1042;runnum<=1042;runnum++){
    beta = new TChain("beta");
    beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_tofcor/correlation2bestT/beta%i*.root",runnum));
    Beta *fbeta = new Beta;
    beta->SetBranchAddress("Beta",&fbeta);
    long n = beta->GetEntries();
    if(n==0) continue;
    long x=0;

    int c0 = 0; // # of imp = 31/30Ne
    int c1 = 0; // # of imp(31/30Ne) with bad GetPixel()
    int c2 = 0; 
    int c3 = 0;
    int c4 = 0;
    int c5 = 0;

    for(x=0;x<n;x++){
      beta->GetEntry(x);
      Implant fimp = fbeta->fImplant;
      if(Ne->IsInside(fimp.fI2S, fimp.fPIN1E)){
        if((fimp.fI2S>11100 && fimp.fI2S<13000) && (fimp.fI2S_I2N>8000 && fimp.fI2S_I2N<16000)){
          c0++;
          pixel piximp = fimp.GetPixel();
          if(fimp.GetPixel()!=piximp) c1++;
          
          for(int m=0;m<fbeta->DecaySize();m++){
            bool flag = false;
            for(auto &it:fbeta->fDecay[m].fGe){
              if(it.GetEnergy()>10 && it.GetEnergy()<4000){
                falg = true;
                break;
              }
            }
            if(!flag) continue;
            c2++;
            pixel pixdec = fbeta->fDecay[m].GetPixel();
            if(fbeta->fDecay[m].GetPixel()!=pixdec) c3++;
          }
        }
      } 



      if((x%5000)==0){
        printf("  on entry %lu / %lu  \r", x,n);
        fflush(stdout);
      }

    }
  }
  printf("on entry %lu / %lu \n", x,n);
  ofile.close(); 

  return;
}
