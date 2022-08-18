



/*
TFile *f = TFile::Open("/home/zhu/packages/BCSSort/gates_add_sub.root");
TFile *f1 = TFile::Open("/home/zhu/packages/BCSSort/gates_add_sub_ab.root");
TFile *f2 = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");

TH1D *single = (TH1D *)f->Get("single");
TH1D *px = (TH1D *)f->Get("px");
TH1D *gadd1 = (TH1D *)f->Get("gadd1");
TH2D *ggmat = (TH2D *)f->Get("ggmat");
TH1D *tp = ggmat->ProjectionX("total_proj",0,4000);

TH1D *gadd1_ab = (TH1D *)f1->Get("gadd1_ab");
TH2D *ggmat_ab = (TH2D *)f1->Get("ggmat_ab");
TH1D *tp_ab = ggmat_ab->ProjectionX("total_proj_ab",0,4000);

TH2D *ggmat_ab_70_32Na = (TH2D *)f2->Get("ggmat_ab_noct_70_32Na");
TH1D *ab_px = ggmat_ab_70_32Na->ProjectionX("ggmat_ab_noct_70ms_32Na_px", 0,4000);
*/

void Process(TH1 *hist){

  std::ofstream ofile;
  
  ofile.open(Form("%s.dat",hist->GetName()));

  long n = hist->GetNbinsX(); 
  long x = 1; // bin = 0 ==> underflow

  for(x=1;x<=n;x++){
    double binwidth   = hist->GetBinWidth(x);
    double bincenter  = hist->GetBinCenter(x);
    double bincontent = hist->GetBinContent(x);
      
    ofile << bincenter-binwidth/2 << "\t" << bincontent << std::endl;  
  }

  ofile.close();
}


void convert_TH1_txt(TH1 *px){
  px->Rebin(1);
  //Process(ab_px);
  Process(px);
}
