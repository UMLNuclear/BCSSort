

{

  TGraphErrors *gr[10];
  for(int i=0;i<10;i++){
    gr[i] = new TGraphErrors(Form("/home/zhu/packages/BCSSort/CORAnalysis/373_fitdat/Cut_31Ne/bin%ims_30ms.dat",i+6));
    gr[i]->SetName(Form("bin%ims",i+6));
    gr[i]->SetTitle(Form("31Ne Addback bin in %.1fms",(i+6)/10.));
  }

  double tlimit = 30.0; // time window;
  TF1 *fx[10]; // for fitting
  TF1 *fg[10]; // for guessing
  for(int i=0;i<10;i++){
    fx[i] = new TF1(Form("fx%i",i), "[0]*exp(-0.693/[1]*x)+[2]",0,tlimit);
    fx[i]->SetParameters(15,2,2);
    fg[i] = new TF1(Form("fg%i",i), "[0]*exp(-0.693/[1]*x)+[2]",0,tlimit);
    fg[i]->SetLineColor(kBlue);
    fg[i]->SetLineWidth(2);
  }
  
  TCanvas *c = new TCanvas;
  c->Divide(4,3);
  for(int i=0;i<10;i++){
    gr[i]->GetXaxis()->SetRangeUser(0,tlimit);
    c->cd(i+1);
    gr[i]->Draw("A*");
    printf("\n bin in %i ms\n",i+6);
    gr[i]->Fit(fx[i]);
    fg[i]->SetParameters(fx[i]->GetParameters());
    fg[i]->Draw("same");
    TText *p0 = new TText;
    p0->SetNDC();
    p0->SetTextFont(1);
    p0->SetTextColor(1);
    p0->SetTextSize(0.05);
    p0->DrawText(0.5,0.8,Form("p0 = %.5f",fx[i]->GetParameter(0)));
    TText *p1 = new TText;
    p1->SetNDC();
    p1->SetTextFont(1);
    p1->SetTextColor(1);
    p1->SetTextSize(0.05);
    p1->DrawText(0.5,0.7,Form("p1 = %.5f",fx[i]->GetParameter(1)));
    TText *p2 = new TText;
    p2->SetNDC();
    p2->SetTextFont(1);
    p2->SetTextColor(1);
    p2->SetTextSize(0.05);
    p2->DrawText(0.5,0.6,Form("p2 = %.5f",fx[i]->GetParameter(2)));
    c->Update();
  } 

  for(int i=0;i<10;i++){
  }

}
