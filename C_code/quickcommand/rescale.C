{
  TH2D *dt_singles_ne31l = (TH2D *)_file0->Get("dt100_addback_ne31l");
  TH2D *dt_singles_ne31r = (TH2D *)_file0->Get("dt100_addback_ne31r");
  TH2D *dt_singles_ne30 = (TH2D *)_file1->Get("dt100_addback_ne30");
  TH2D *dt_singles_ne31 = (TH2D *)_file2->Get("dt100_addback_ne31");

  TH1D *h0  = (TH1D *)dt_singles_ne30->ProjectionY("ne30_10ms"  ,1,10);
  TH1D *h1  = (TH1D *)dt_singles_ne31->ProjectionY("ne31_10ms"  ,1,10);
  TH1D *h1l = (TH1D *)dt_singles_ne31l->ProjectionY("ne31l_10ms",1,10);
  TH1D *h1r = (TH1D *)dt_singles_ne31r->ProjectionY("ne31r_10ms",1,10);

  TH1D *ch0  = (TH1D *)h0->Clone(Form("%s_clone",h0->GetName()));
  TH1D *ch1  = (TH1D *)h1->Clone(Form("%s_clone",h1->GetName()));
  TH1D *ch1l = (TH1D *)h1l->Clone(Form("%s_clone",h1l->GetName()));
  TH1D *ch1r = (TH1D *)h1r->Clone(Form("%s_clone",h1r->GetName()));

  int energy[5] = {150, 365, 1597, 1963, 2114};
  //intscl = integral for each peak shown in energy[5] from addback within 10ms.(FAILED)
  double intscl[5][3] = {{0.284789, 0.051355, 0.233133},
                         {0.242547, 0.048780, 0.196477},
                         {0.228426, 0.025381, 0.203046},
                         {0.286713, 0.048951, 0.244755},
                         {0.305882, 0.070588, 0.235294}};


  int numE = 1;

  ch1 ->Add(h0,-intscl[numE][0]);
  ch1l->Add(h0,-intscl[numE][1]);
  ch1r->Add(h0,-intscl[numE][2]);

  ch0 ->Rebin(2);
  ch1 ->Rebin(2);
  ch1l->Rebin(2);
  ch1r->Rebin(2);

  ch0 -> SetTitle(Form("Ne30 Addback 10ms rescale@%i",energy[numE]));
  ch1 -> SetTitle(Form("Ne31 Addback 10ms rescale@%i",energy[numE]));
  ch1l-> SetTitle(Form("Ne31L Addback 10ms rescale@%i",energy[numE]));
  ch1r-> SetTitle(Form("Ne31R Addback 10ms rescale@%i",energy[numE]));

  TCanvas *c = new TCanvas;
  c->Divide(1,4);
  c->cd(1); ch0 ->Draw();
  c->cd(2); ch1 ->Draw();
  c->cd(3); ch1l->Draw();
  c->cd(4); ch1r->Draw();
  ch0 -> SetLineColor(kRed);
  ch1 -> SetLineColor(kRed);
  ch1l-> SetLineColor(kRed);
  ch1r-> SetLineColor(kRed);

  h0 ->Rebin(2);
  h1 ->Rebin(2);
  h1l->Rebin(2);
  h1r->Rebin(2);
  TCanvas *c1 = new TCanvas;
  c1->Divide(1,4);
  c1->cd(1); h0 ->Draw(); ch0 ->Draw("same");
  c1->cd(2); h1 ->Draw(); ch1 ->Draw("same");
  c1->cd(3); h1l->Draw(); ch1l->Draw("same");
  c1->cd(4); h1r->Draw(); ch1r->Draw("same");
}
