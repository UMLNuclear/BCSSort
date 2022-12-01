





TH1D *thing(double x){

  TH1D *ch1 = (TH1D*)h1->Clone("ch1");  //31Na
  TH1D *ch0 = (TH1D*)h0->Clone("ch0");  //30Na

  ch1->GetXaxis()->SetRangeUser(0,4000);
  ch0->GetXaxis()->SetRangeUser(0,4000);
  TH1 *bg1 = TSpectrum::StaticBackground(ch1);
  TH1 *bg0 = TSpectrum::StaticBackground(ch0);

  ch1->Add(bg1,-1);
  ch0->Add(bg0,-1);
 
  //x = 0.365; // for 10ms 
  x = 0.345;  // for 100ms
    

  TH1D *cch1 = (TH1D *)ch1->Clone("cch1");
  cch1->Add(ch0,-x);
  ch1->Draw();  
  ch0->Scale(x);
  ch0->SetLineColor(kRed);  
  ch0->SetLineWidth(2);  
  ch0->Sumw2(0); 
  ch0->Draw("same hist");
 
  return cch1;


}
