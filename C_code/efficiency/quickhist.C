/*{

  //.L C_code/efficiency/efficiency.C
  TH1D *sum = (TH1D *)single_Cal_48->Clone("sum0");
  sum->Add(single_Cal_49,1);
  sum->Add(single_Cal_50,1);
  sum->Add(single_Cal_51,1);
  TSpectrum s1;
  TH1D *csum = (TH1D *)sum->Clone("sum0_clone");
  sbg = (TH1D *)s1.Background(sum,10);
  csum->Add(sbg,-1);
  //nsum = (TH1D *)NegRemove(csum);

  TH1D *adb = (TH1D *)Addback_1->Clone("addback1");
  TH1D *cadb = (TH1D *)Addback_1->Clone("addback1_clone");
  TSpectrum s2;
  abg = (TH1D *)s2.Background(adb,10);
  cadb->Add(abg,-1);
  //nabd = (TH1D *)NegRemove(cabd);

}*/

TH1D *xtal(int xtal=0){

  TH1D *sp = (TH1D *)_file0->Get(Form("single_Cal_%i",xtal));
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_clone",sp->GetName()));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  csp->Add(bg,-1);
  return csp;

}

TH1D *addback(int xtal=0){

  TH1D *sp = (TH1D *)_file0->Get(Form("Addback_%i",xtal));
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_clone",sp->GetName()));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  csp->Add(bg,-1);
  return csp;

}

TH1D *sum(int xtal=0){

  xtal = xtal*4;
  TH1D *temp = (TH1D *)_file0->Get(Form("single_Cal_%i",xtal));
  TH1D *sp = (TH1D *)temp->Clone(Form("sum%i",xtal));
  for(int i=xtal+1;i<xtal+4;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)_file0->Get(Form("single_Cal_%i",i));
    sp->Add(temp,1);  
  }
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_clone",sp->GetName()));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  csp->Add(bg,-1);
  return csp;
}

TH1D *allsum(){

  TH1D *sp = (TH1D *)single_Cal_0->Clone("all_sum");
  for(int i=1;i<64;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)_file0->Get(Form("single_Cal_%i",i));
    sp->Add(temp,1);  
  }
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_clone",sp->GetName()));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  csp->Add(bg,-1);
  return csp;
}

TH1D *alladb(){

  TH1D *sp = (TH1D *)Addback_0->Clone("all_addback");
  for(int i=1;i<16;i++){
    TH1D *temp = (TH1D *)_file0->Get(Form("Addback_%i",i));
    sp->Add(temp,1);  
  }
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_clone",sp->GetName()));
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //csp->Draw();
  //bg->SetLineColor(kRed);
  //bg->Draw("same");
  csp->Add(bg,-1);
  return csp;
}
