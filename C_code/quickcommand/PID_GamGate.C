{

  int energy[7] = {150, 365, 373, 410, 1597, 1963, 2114};
  TH2D *pids10[7];
  TH2D *pids10bg[7];
  TH2D *pids20[7];
  TH2D *pids20bg[7];
  TH2D *pida10[7];
  TH2D *pida10bg[7];
  TH2D *pida20[7];
  TH2D *pida20bg[7];
 
  TH1D *s10[7];
  TH1D *s10bg[7];
  TH1D *s20[7];
  TH1D *s20bg[7];
  TH1D *a10[7];
  TH1D *a10bg[7];
  TH1D *a20[7];
  TH1D *a20bg[7];
 
  TH1D *cs10[7];
  TH1D *cs20[7];
  TH1D *ca10[7];
  TH1D *ca20[7];
  
  TCutG *cut0 = new TCutG("cut0",5);
  cut0->SetPoint(0, 10800, 6400);
  cut0->SetPoint(1, 10800, 4800);
  cut0->SetPoint(2, 15500, 3800);
  cut0->SetPoint(3, 15500, 5400);
  cut0->SetPoint(4, 10800, 6400);
  cut0->SetLineColor(kRed);
  cut0->SetLineWidth(2);

  TH2D *pid = (TH2D *)_file0->Get("pid");
  pid->GetXaxis()->SetRangeUser(8800,15800);
  pid->GetYaxis()->SetRangeUser(4800,7800);
  new TCanvas;
  pid->Draw("colz");
  cut0->Draw("same");
  
  TLine *l1 = new TLine(10970,0,10970,700); // left  boundary or 31Ne
  TLine *l2 = new TLine(11470,0,11470,700); // Mid line of 31Ne
  TLine *l3 = new TLine(11970,0,11970,700); // left  boundary of 30Ne
  TLine *l4 = new TLine(12870,0,12870,700); // right boundary of 30Ne

  l1->SetLineColor(kRed);
  l2->SetLineColor(kRed);
  l3->SetLineColor(kRed);
  l4->SetLineColor(kRed);
  l1->SetLineWidth(3);
  l2->SetLineWidth(3);
  l3->SetLineWidth(3);
  l4->SetLineWidth(3);


  TCanvas *c[7];
  for(int i=0;i<7;i++){
    pids10[i]   = (TH2D *)_file0->Get(Form("pid_singles_10ms_%i",energy[i]));
    pids10bg[i] = (TH2D *)_file0->Get(Form("pid_singles_10ms_%iBG",energy[i]));
    pids20[i]   = (TH2D *)_file0->Get(Form("pid_singles_20ms_%i",energy[i]));
    pids20bg[i] = (TH2D *)_file0->Get(Form("pid_singles_20ms_%iBG",energy[i]));
    pida10[i]   = (TH2D *)_file0->Get(Form("pid_addback_10ms_%i",energy[i]));
    pida10bg[i] = (TH2D *)_file0->Get(Form("pid_addback_10ms_%iBG",energy[i]));
    pida20[i]   = (TH2D *)_file0->Get(Form("pid_addback_20ms_%i",energy[i]));
    pida20bg[i] = (TH2D *)_file0->Get(Form("pid_addback_20ms_%iBG",energy[i]));
  
    s10[i]   = pids10[i]  ->ProjectionX(Form("tof_singles_10ms_%i",energy[i]),  1,720,"[cut0]");
    s10bg[i] = pids10bg[i]->ProjectionX(Form("tof_singles_10ms_%iBG",energy[i]),1,720,"[cut0]");
    s20[i]   = pids20[i]  ->ProjectionX(Form("tof_singles_20ms_%i",energy[i]),  1,720,"[cut0]");
    s20bg[i] = pids20bg[i]->ProjectionX(Form("tof_singles_20ms_%iBG",energy[i]),1,720,"[cut0]");
    a10[i]   = pida10[i]  ->ProjectionX(Form("tof_addback_10ms_%i",energy[i]),  1,720,"[cut0]");
    a10bg[i] = pida10bg[i]->ProjectionX(Form("tof_addback_10ms_%iBG",energy[i]),1,720,"[cut0]");
    a20[i]   = pida20[i]  ->ProjectionX(Form("tof_addback_20ms_%i",energy[i]),  1,720,"[cut0]");
    a20bg[i] = pida20bg[i]->ProjectionX(Form("tof_addback_20ms_%iBG",energy[i]),1,720,"[cut0]");
   
    s10[i]->SetTitle(Form("PID singles 10ms gated %i",energy[i]));
    s20[i]->SetTitle(Form("PID singles 20ms gated %i",energy[i]));
    a10[i]->SetTitle(Form("PID addback 10ms gated %i",energy[i]));
    a20[i]->SetTitle(Form("PID addback 20ms gated %i",energy[i]));
 
    cs10[i] = (TH1D *)s10[i]->Clone(Form("%s_cl",s10[i]->GetName()));
    cs20[i] = (TH1D *)s20[i]->Clone(Form("%s_cl",s20[i]->GetName()));
    ca10[i] = (TH1D *)a10[i]->Clone(Form("%s_cl",a10[i]->GetName()));
    ca20[i] = (TH1D *)a20[i]->Clone(Form("%s_cl",a20[i]->GetName()));

    cs10[i]->Add(s10bg[i],-1);   
    cs20[i]->Add(s20bg[i],-1);   
    ca10[i]->Add(a10bg[i],-1);   
    ca20[i]->Add(a20bg[i],-1);   
    cs10[i]->SetLineColor(kBlue);   
    cs20[i]->SetLineColor(kBlue);   
    ca10[i]->SetLineColor(kBlue);   
    ca20[i]->SetLineColor(kBlue);   
    cs10[i]->SetLineWidth(2);   
    cs20[i]->SetLineWidth(2);   
    ca10[i]->SetLineWidth(2);   
    ca20[i]->SetLineWidth(2);   

    s10[i]->Rebin(4);
    s20[i]->Rebin(4);
    a10[i]->Rebin(4);
    a20[i]->Rebin(4);
    cs10[i]->Rebin(4);
    cs20[i]->Rebin(4);
    ca10[i]->Rebin(4);
    ca20[i]->Rebin(4);

    c[i] = new TCanvas;
    c[i]->Divide(2,3);
    c[i]->cd(1);
    pida10[i]->Draw("colz");
    cut0->Draw("same");
    c[i]->cd(2);
    pida20[i]->Draw("colz");
    cut0->Draw("same");
    c[i]->cd(3);
    s10[i]->Draw();
    cs10[i]->Draw("same");
    l1->Draw("same");l2->Draw("same");l3->Draw("same");l4->Draw("same");
    c[i]->cd(4);
    s20[i]->Draw();
    cs20[i]->Draw("same");
    l1->Draw("same");l2->Draw("same");l3->Draw("same");l4->Draw("same");
    c[i]->cd(5);
    a10[i]->Draw();
    ca10[i]->Draw("same");
    l1->Draw("same");l2->Draw("same");l3->Draw("same");l4->Draw("same");
    c[i]->cd(6);
    a20[i]->Draw();
    ca20[i]->Draw("same");
    l1->Draw("same");l2->Draw("same");l3->Draw("same");l4->Draw("same");
    c[i]->Update();  
    
  }


}
