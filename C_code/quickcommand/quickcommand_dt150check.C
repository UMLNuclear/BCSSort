{
  TH2D *dt_singles_ne31l = (TH2D *)_file0->Get("dt30_singles_ne31l");
  TH2D *dt_singles_ne31r = (TH2D *)_file0->Get("dt30_singles_ne31r");
  TH2D *dt_singles_ne30  = (TH2D *)_file1->Get("dt30_singles_ne30");
  TH2D *dt_singles_ne31  = (TH2D *)_file2->Get("dt30_singles_ne31");

  TH1D *h0  = (TH1D *)dt_singles_ne30 ->ProjectionX("ne30", 147,152);
  TH1D *h1  = (TH1D *)dt_singles_ne31 ->ProjectionX("ne31", 147,152);
  TH1D *h1l = (TH1D *)dt_singles_ne31l->ProjectionX("ne31l",147,152);
  TH1D *h1r = (TH1D *)dt_singles_ne31r->ProjectionX("ne31r",147,152);
  TH1D *h[4];
  h[0] = (TH1D *)h0 -> Clone(Form("%s_c",h0 ->GetName()));
  h[1] = (TH1D *)h1 -> Clone(Form("%s_c",h1 ->GetName()));
  h[2] = (TH1D *)h1l-> Clone(Form("%s_c",h1l->GetName()));
  h[3] = (TH1D *)h1r-> Clone(Form("%s_c",h1r->GetName()));

  TH1D *bg0  = (TH1D *)dt_singles_ne30 ->ProjectionX("bgne30", 154,184);
  TH1D *bg1  = (TH1D *)dt_singles_ne31 ->ProjectionX("bgne31", 154,184);
  TH1D *bg1l = (TH1D *)dt_singles_ne31l->ProjectionX("bgne31l",154,184);
  TH1D *bg1r = (TH1D *)dt_singles_ne31r->ProjectionX("bgne31r",154,184);
  TH1D *bg[4];
  bg[0] = (TH1D *)bg0 -> Clone(Form("%s_c",bg0 ->GetName()));
  bg[1] = (TH1D *)bg1 -> Clone(Form("%s_c",bg1 ->GetName()));
  bg[2] = (TH1D *)bg1l-> Clone(Form("%s_c",bg1l->GetName()));
  bg[3] = (TH1D *)bg1r-> Clone(Form("%s_c",bg1r->GetName()));
  
  TF1 *fx[4];
  TCanvas *c = new TCanvas;
  c->Divide(3,1);
  for(int i=1;i<4;i++){
    bg[i]->Scale(1/6.);
    h[i]->Add(bg[i],-1);
    h[i]->Rebin(5);
    h[i]->GetXaxis()->SetRangeUser(0,30);
    h[i]->SetTitle(Form("decaytime @150keV %s",h[i]->GetName()));
    c->cd(i);
    h[i]->Draw("hist");

    fx[i] = new TF1(Form("fx%i",i), "[0]*exp(-0.69314718/[1]*x)+[2]",0,30);
    fx[i]->SetParameters(50,4,0.1);
    fx[i]->SetParLimits(0,0,200);
    fx[i]->SetParLimits(1,0,10);
    fx[i]->SetParLimits(2,0,50);

    h[i]->Fit(fx[i]);
    TF1 *fit = new TF1("fit","[0]*exp(-0.69314718/[1]*x)+[2]",0,30);
    fit->SetParameters(fx[i]->GetParameter(0),fx[i]->GetParameter(1),fx[i]->GetParameter(2));
    fit->SetLineColor(kRed);
    fit->Draw("same");
    
    c->Update();
    double y = gPad->GetUymax();
    TLatex l;
    l.SetTextSize(0.05);
    l.DrawLatex(13.5,y*.45,Form("t_{1/2} = %.3f(%.3f)ms",fx[i]->GetParameter(1), fx[i]->GetParError(1)));
 
  }
    

}
