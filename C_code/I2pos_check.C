

#include<util.h>



void I2pos_check(){


  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");
  long n0=0;
  long n=0;
  long x=0;
  for(int runnum=1023;runnum<1140;runnum++){
    TChain *beta = new TChain("beta");
    beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_exdTOF_tofcor/correlation1bestT/beta%i*.root",runnum));
    n = beta->GetEntries();
    if(n==0) continue;
    Beta *fbeta = new Beta;
    beta->SetBranchAddress("Beta",&fbeta);
  
    for(x=0;x<n;x++){
      beta->GetEntry(x);
      Implant fimp = fbeta->fImplant;
      if(Ne->IsInside(fimp.fI2S, fimp.fPIN1E)){
        FillHistogram("I2pos",6000,0,6000000,x+n0,800,0,40000,fimp.fI2S_I2N);
      }
    }
    n0 += n;
  }
  SaveHistograms("I2pos_check.root");
  return;
}
