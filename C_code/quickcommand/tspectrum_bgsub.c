//=============================================//
//3 files + addback + 10ms;
//Spectrum + TSpectrum BG subtraction
//=============================================//
{
  TH2D *dt_singles_ne31l = (TH2D *)_file0->Get("dt100_addback_ne31l");
  TH2D *dt_singles_ne31r = (TH2D *)_file0->Get("dt100_addback_ne31r");
  TH2D *dt_singles_ne30 = (TH2D *)_file1->Get("dt100_addback_ne30");
  TH2D *dt_singles_ne31 = (TH2D *)_file2->Get("dt100_addback_ne31");

  TH1D *h0  = (TH1D *)dt_singles_ne30->ProjectionY("ne30_10ms",1,10);
  TH1D *h1  = (TH1D *)dt_singles_ne31->ProjectionY("ne31_10ms",1,10);
  TH1D *h1l = (TH1D *)dt_singles_ne31l->ProjectionY("ne31l_10ms",1,10);
  TH1D *h1r = (TH1D *)dt_singles_ne31r->ProjectionY("ne31r_10ms",1,10);

  TH1D *ch0  = (TH1D *)h0->Clone(Form("%s_clone",h0->GetName()));
  TH1D *ch1  = (TH1D *)h1->Clone(Form("%s_clone",h1->GetName()));
  TH1D *ch1l = (TH1D *)h1l->Clone(Form("%s_clone",h1l->GetName()));
  TH1D *ch1r = (TH1D *)h1r->Clone(Form("%s_clone",h1r->GetName()));

  TSpectrum s;
  TH1D *bg[4];
  bg[0] = (TH1D *)s.Background(h0,25);
  bg[1] = (TH1D *)s.Background(h1,25);
  bg[2] = (TH1D *)s.Background(h1l,25);
  bg[3] = (TH1D *)s.Background(h1r,25);

  ch0 ->Add(bg[0],-1);
  ch1 ->Add(bg[1],-1);
  ch1l->Add(bg[2],-1);
  ch1r->Add(bg[3],-1);

  //ch0 ->Rebin(2);
  //ch1 ->Rebin(2);
  //ch1l->Rebin(2);
  //ch1r->Rebin(2);

  TCanvas *c1 = new TCanvas;
  c1->Divide(1,4);
  c1->cd(1); h0 ->Draw();
  c1->cd(2); h1 ->Draw();
  c1->cd(3); h1l->Draw();
  c1->cd(4); h1r->Draw();
  
  TCanvas *c = new TCanvas;
  c->Divide(1,4);
  c->cd(1); ch0 ->Draw();
  c->cd(2); ch1 ->Draw();
  c->cd(3); ch1l->Draw();
  c->cd(4); ch1r->Draw();
}
