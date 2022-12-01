void quickcommand(TH1 *h){

  TF1 *fx  = new TF1("fx", "[0]*exp(-0.693/[1]*x)+[2]*exp(-0.693/[3]*x)+[4]");
  TF1 *fx1 = new TF1("fx1","[0]*exp(-0.693/[1]*x)+[2]");
  TF1 *fx2 = new TF1("fx2","[0]*exp(-0.693/[1]*x)+[2]");

  fx->SetParameters(40,7.9,10,2,0);
  fx->SetParLimits(0,0,100000);
  fx->SetParLimits(1,5,10);
  fx->SetParLimits(2,0,1000);
  fx->SetParLimits(3,0,5);
  fx->SetParLimits(4,0,100000);
  fx->FixParameter(1,7.9);
  fx->FixParameter(3,1.1);

  h->Fit(fx);
  fx1->SetParameter(0,fx->GetParameter(0));
  fx1->SetParameter(1,fx->GetParameter(1));
  fx1->SetParameter(2,fx->GetParameter(4));
  fx2->SetParameter(0,fx->GetParameter(2));
  fx2->SetParameter(1,fx->GetParameter(3));
  fx2->SetParameter(2,fx->GetParameter(4));

  fx1->SetLineColor(kBlue);
  fx2->SetLineColor(kGreen);
  fx1->SetRange(0,40);
  fx2->SetRange(0,40);
}
