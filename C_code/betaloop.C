#include<util.h>
#include<globals.h>

void betaloop(){

  TChannel::ReadDetMapFile();

  long c0 = 0;
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

  std::map<pixel, double>  fImpTMap;
  std::map<pixel, Implant*>fImpMap;

  int dnum;
  double dt, dt_nb;

  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_Na32cut.root");
  TCutG *Na32 = (TCutG *)cutf->Get("Na32");
  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *Ne = (TCutG *)cutf1->Get("Ne");

  for(int runnum=1040;runnum<=1040;runnum++){
    TChain *beta = new TChain("beta");
    beta->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/beta/beta_good_prompt_tofcor/beta%i*.root",runnum));
    long n = beta->GetEntries();
    if(n==0) continue;
    c5+=n;
    Beta *fbeta = new Beta;
    beta->SetBranchAddress("Beta",&fbeta);
    
    long x =0;

    for(int i=0;i<40;i++){
      for(int j=0;j<40;j++){
        pixel pix = std::make_pair(i,j);
        fImpTMap[pix] = -1;
        fImpMap[pix]  = new Implant;
      }
    }

    //n = 1e2;
    for(x=0;x<n;x++){
      beta->GetEntry(x);
      Implant fimp = fbeta->fImplant;
      if(Ne->IsInside(fimp.fI2S, fimp.fPIN1E)){
        if(fimp.fI2S_I2N<8000 || fimp.fI2S_I2N>16000) continue;
        FillHistogram("decaysize_50ms", 100,0,100,fbeta->DecaySize(50));
      }
    }

  //std::cout<<"total entries(implant) = " << c5 << std::endl;
  //std::cout<<"implant(decaysize==0) = " << c0 <<std::endl; 
  //std::cout<<"implant(SSSD>0) = "       << c1 <<std::endl; 
  //std::cout<<std::endl;
  //std::cout<<"32Na imp               = "<< c2 <<std::endl; 
  //std::cout<<"32Na imp(decaysize==0) = "<< c3 <<std::endl; 
  //std::cout<<"32Na imp(SSSD>0) = "      << c4 <<std::endl; 

  SaveHistograms("betaloop_op.root");

  return;
}
