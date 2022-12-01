



void histsort(){


  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/pid_cut.root");
  TCutG *na = (TCutG *)cutf->Get("Na");
  TCutG *ne = (TCutG *)cutf->Get("Ne");
  na->SetName("na");
  ne->SetName("ne");

  TFile *hist1 = TFile::Open("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/cor1_hist.root");
  TFile *hist2 = TFile::Open("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/cor2_hist.root");
  TH2D *p1[6];
  TH2D *b1[6];
  TH2D *p2[6];
  TH2D *b2[6];
  
  TH1D *nap1[6];
  TH1D *nab1[6];
  TH1D *nap2[6];
  TH1D *nab2[6];
  TH1D *nep1[6];
  TH1D *neb1[6];
  TH1D *nep2[6];
  TH1D *neb2[6];
  
  TH1D *nac1[6];
  TH1D *nac2[6];
  TH1D *nec1[6];
  TH1D *nec2[6];

  for(int i=0;i<6;i++){
    p1[i] = (TH2D *)hist1->Get(Form("pidg%i",i)); 
    b1[i] = (TH2D *)hist1->Get(Form("pidb%i",i)); 
    p2[i] = (TH2D *)hist2->Get(Form("pidg%i",i)); 
    b2[i] = (TH2D *)hist2->Get(Form("pidb%i",i)); 
  }

  for(int i=0;i<6;i++){
    nap1[i] = p1[i]->ProjectionX(Form("%s_na",p1[i]->GetName()),1,300,"[na]");
    nab1[i] = b1[i]->ProjectionX(Form("%s_na",b1[i]->GetName()),1,300,"[na]");
    nap2[i] = p2[i]->ProjectionX(Form("%s_na",p2[i]->GetName()),1,300,"[na]");
    nab2[i] = b2[i]->ProjectionX(Form("%s_na",b2[i]->GetName()),1,300,"[na]");
    nep1[i] = p1[i]->ProjectionX(Form("%s_ne",p1[i]->GetName()),1,300,"[ne]");
    neb1[i] = b1[i]->ProjectionX(Form("%s_ne",b1[i]->GetName()),1,300,"[ne]");
    nep2[i] = p2[i]->ProjectionX(Form("%s_ne",p2[i]->GetName()),1,300,"[ne]");
    neb2[i] = b2[i]->ProjectionX(Form("%s_ne",b2[i]->GetName()),1,300,"[ne]");

    nac1[i] = (TH1D *)nap1[i]->Clone(Form("c%s",nap1[i]->GetName()));
    nac2[i] = (TH1D *)nap2[i]->Clone(Form("c%s",nap2[i]->GetName()));
    nec1[i] = (TH1D *)nep1[i]->Clone(Form("c%s",nep1[i]->GetName()));
    nec2[i] = (TH1D *)nep2[i]->Clone(Form("c%s",nep2[i]->GetName()));

    nac1[i]->Add(nab1[i], -1);
    nac2[i]->Add(nab2[i], -1);
    nec1[i]->Add(neb1[i], -1);
    nec2[i]->Add(neb2[i], -1);

  }

  TFile *newf = new TFile("/home/zhu/packages/BCSSort/CORAnalysis/example/hist/histsort_corl12.root","recreate");
  for(int i=0;i<6;i++){
    nap1[i] -> Write();
    nab1[i] -> Write();
    nap2[i] -> Write();
    nab2[i] -> Write();
    nep1[i] -> Write();
    neb1[i] -> Write();
    nep2[i] -> Write();
    neb2[i] -> Write();
    nac1[i] -> Write();
    nac2[i] -> Write();
    nec1[i] -> Write();
    nec2[i] -> Write();
  }
  newf->Close();

}










