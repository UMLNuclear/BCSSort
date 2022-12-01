


{
  
  TH2D *h1 = (TH2D *)_file0->Get("pid1_1");
  TH2D *h2 = (TH2D *)_file0->Get("pid4_1");
  TH2D *h3 = (TH2D *)_file0->Get("pid6_1");
  TH2D *h4 = (TH2D *)_file0->Get("pid8_1");
  TH2D *h5 =(TH2D *)_file0->Get("pid10_1");
  TH2D *b1 = (TH2D *)_file0->Get("pid3_1");
  TH2D *b2 = (TH2D *)_file0->Get("pid5_1");
  TH2D *b3 = (TH2D *)_file0->Get("pid7_1");
  TH2D *b4 = (TH2D *)_file0->Get("pid9_1");
  TH2D *b5 =(TH2D *)_file0->Get("pid11_1");

  TFile *cutf = new TFile("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *na = (TCutG *)cutf->Get("Na");
  TCutG *ne = (TCutG *)cutf->Get("Ne");
  na->SetName("na");
  ne->SetName("ne");

  TH1D *h1na = (TH1D *)h1->ProjectionX("h1na",1,300,"[na]");
  TH1D *h2na = (TH1D *)h2->ProjectionX("h2na",1,300,"[na]");
  TH1D *h3na = (TH1D *)h3->ProjectionX("h3na",1,300,"[na]");
  TH1D *h4na = (TH1D *)h4->ProjectionX("h4na",1,300,"[na]");
  TH1D *h5na = (TH1D *)h5->ProjectionX("h5na",1,300,"[na]");
  TH1D *h1ne = (TH1D *)h1->ProjectionX("h1ne",1,300,"[ne]");
  TH1D *h2ne = (TH1D *)h2->ProjectionX("h2ne",1,300,"[ne]");
  TH1D *h3ne = (TH1D *)h3->ProjectionX("h3ne",1,300,"[ne]");
  TH1D *h4ne = (TH1D *)h4->ProjectionX("h4ne",1,300,"[ne]");
  TH1D *h5ne = (TH1D *)h5->ProjectionX("h5ne",1,300,"[ne]");
  TH1D *b1na = (TH1D *)b1->ProjectionX("bg1na",1,300,"[na]");
  TH1D *b2na = (TH1D *)b2->ProjectionX("bg2na",1,300,"[na]");
  TH1D *b3na = (TH1D *)b3->ProjectionX("bg3na",1,300,"[na]");
  TH1D *b4na = (TH1D *)b4->ProjectionX("bg4na",1,300,"[na]");
  TH1D *b5na = (TH1D *)b5->ProjectionX("bg5na",1,300,"[na]");
  TH1D *b1ne = (TH1D *)b1->ProjectionX("bg1ne",1,300,"[ne]");
  TH1D *b2ne = (TH1D *)b2->ProjectionX("bg2ne",1,300,"[ne]");
  TH1D *b3ne = (TH1D *)b3->ProjectionX("bg3ne",1,300,"[ne]");
  TH1D *b4ne = (TH1D *)b4->ProjectionX("bg4ne",1,300,"[ne]");
  TH1D *b5ne = (TH1D *)b5->ProjectionX("bg5ne",1,300,"[ne]");
 
  h1na->Rebin(8);
  h2na->Rebin(8);
  h3na->Rebin(8);
  h4na->Rebin(8);
  h5na->Rebin(8);
  h1ne->Rebin(8);
  h2ne->Rebin(8);
  h3ne->Rebin(8);
  h4ne->Rebin(8);
  h5ne->Rebin(8);
  b1na->Rebin(8);
  b2na->Rebin(8);
  b3na->Rebin(8);
  b4na->Rebin(8);
  b5na->Rebin(8);
  b1ne->Rebin(8);
  b2ne->Rebin(8);
  b3ne->Rebin(8);
  b4ne->Rebin(8);
  b5ne->Rebin(8);

 
  TH1D *ch1na = (TH1D *)h1na->Clone("tof1na"); 
  TH1D *ch2na = (TH1D *)h2na->Clone("tof2na"); 
  TH1D *ch3na = (TH1D *)h3na->Clone("tof3na"); 
  TH1D *ch4na = (TH1D *)h3na->Clone("tof4na"); 
  TH1D *ch5na = (TH1D *)h3na->Clone("tof5na"); 
  TH1D *ch1ne = (TH1D *)h1ne->Clone("tof1ne"); 
  TH1D *ch2ne = (TH1D *)h2ne->Clone("tof2ne"); 
  TH1D *ch3ne = (TH1D *)h3ne->Clone("tof3ne"); 
  TH1D *ch4ne = (TH1D *)h3na->Clone("tof4ne"); 
  TH1D *ch5ne = (TH1D *)h3na->Clone("tof5ne"); 

  h1na->SetTitle("150keV in Na chain within 10ms"); 
  h2na->SetTitle("373keV in Na chain within 10ms"); 
  h3na->SetTitle("885keV in Na chain within 10ms"); 
  h4na->SetTitle("72keV  in Na chain within 10ms"); 
  h5na->SetTitle("1516keV in Na chain within 10ms"); 
  h1ne->SetTitle("150keV in Ne chain within 10ms"); 
  h2ne->SetTitle("373keV in Ne chain within 10ms"); 
  h3ne->SetTitle("885keV in Ne chain within 10ms"); 
  h4ne->SetTitle("72keV  in Ne chain within 10ms"); 
  h5ne->SetTitle("1516keV in Ne chain within 10ms"); 

  h1na->GetXaxis()->SetRangeUser(10000, 16000); 
  h2na->GetXaxis()->SetRangeUser(10000, 16000); 
  h3na->GetXaxis()->SetRangeUser(10000, 16000); 
  h4na->GetXaxis()->SetRangeUser(10000, 16000); 
  h5na->GetXaxis()->SetRangeUser(10000, 16000); 
  h1ne->GetXaxis()->SetRangeUser(10000, 16000); 
  h2ne->GetXaxis()->SetRangeUser(10000, 16000); 
  h3ne->GetXaxis()->SetRangeUser(10000, 16000); 
  h4ne->GetXaxis()->SetRangeUser(10000, 16000); 
  h5ne->GetXaxis()->SetRangeUser(10000, 16000); 

  ch1na->SetLineColor(kBlue); 
  ch2na->SetLineColor(kBlue); 
  ch3na->SetLineColor(kBlue); 
  ch4na->SetLineColor(kBlue); 
  ch5na->SetLineColor(kBlue); 
  ch1ne->SetLineColor(kBlue); 
  ch2ne->SetLineColor(kBlue); 
  ch3ne->SetLineColor(kBlue); 
  ch4ne->SetLineColor(kBlue); 
  ch5ne->SetLineColor(kBlue); 

  ch1na->SetLineWidth(4); 
  ch2na->SetLineWidth(4); 
  ch3na->SetLineWidth(4); 
  ch1ne->SetLineWidth(4); 
  ch2ne->SetLineWidth(4); 
  ch3ne->SetLineWidth(4); 

  ch1na->Add(b1na,-1); 
  ch2na->Add(b2na,-1); 
  ch3na->Add(b3na,-1); 
  ch4na->Add(b4na,-1); 
  ch5na->Add(b5na,-1); 
  ch1ne->Add(b1ne,-1); 
  ch2ne->Add(b2ne,-1); 
  ch3ne->Add(b3ne,-1); 
  ch4ne->Add(b4ne,-1); 
  ch5ne->Add(b5ne,-1); 

  TCanvas *c1 = new TCanvas;
  c1->Divide(2,1);
  c1->cd(1); h1na->Draw("hist"); ch1na->Draw("same hist");
  c1->cd(2); h1ne->Draw("hist"); ch1ne->Draw("same hist");
  TCanvas *c2 = new TCanvas;
  c2->Divide(2,1);
  c2->cd(1); h2na->Draw("hist"); ch2na->Draw("same hist");
  c2->cd(2); h2ne->Draw("hist"); ch2ne->Draw("same hist");
  TCanvas *c3 = new TCanvas;
  c3->Divide(2,1);
  c3->cd(1); h3na->Draw("hist"); ch3na->Draw("same hist");
  c3->cd(2); h3ne->Draw("hist"); ch3ne->Draw("same hist");
  //TCanvas *c4 = new TCanvas;
  //c4->Divide(2,1);
  //c4->cd(1); h4na->Draw("hist"); ch4na->Draw("same hist");
  //c4->cd(2); h4ne->Draw("hist"); ch4ne->Draw("same hist");
  //TCanvas *c5 = new TCanvas;
  //c5->Divide(2,1);
  //c5->cd(1); h5na->Draw("hist"); ch5na->Draw("same hist");
  //c5->cd(2); h5ne->Draw("hist"); ch5ne->Draw("same hist");
  

}
