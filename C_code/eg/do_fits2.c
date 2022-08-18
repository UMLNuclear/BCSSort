

void do_fits2(int rebin=10) {

  TH2* pw = (TH2*)_file0->Get("dt_TOF_30ms_150");
  //TH2* bg = (TH2*)_file0->Get("dt_TOF_30ms_BGR150");
  TH2* bg = (TH2*)_file0->Get("dt_TOF_30ms_extBGR150");
  bg->Scale(1/6.);
  TH2* gam = (TH2*)_file0->Get("singles30ms_TOF_Ne");

  pw = (TH2*)pw->Clone();
  bg = (TH2*)bg->Clone();

  //pw->RebinX(10);
  //bg->RebinX(10);
  pw->RebinY(rebin);
  bg->RebinY(rebin);

  int low,high;
  double width = 300;
  double shift = 100;
  double current = 10770;

  TH1 *p[20];
  TH1 *b[20];
  TH1 *pb[20];
  TH1 *singles[16];

  TF1 *f[20];

  gStyle->SetOptStat("nei");

  for(int x=0;x<20;x++) {
    low  = pw->GetXaxis()->FindBin(current+1);
    high = pw->GetXaxis()->FindBin(current+width-1);
    //singles[x] = gam->ProjectionY(Form("gam%i",x),low,high);
    p[x] = pw->ProjectionY(Form("p%i",x),low,high); // from 11500 to 12000
    //p[x]->Rebin(10);
    b[x] = bg->ProjectionY(Form("b%i",x),low,high); // from 11500 to 12000
    //b[x]->Rebin(10);
    pb[x] = (TH1*)p[x]->Clone(Form("pb%i",x));
    pb[x]->Add(b[x],-1);
    pb[x]->SetTitle(Form("%i to %i",(int)current,int(current+width)));
    current+=shift;
  }





/*
  TF1 *f1  = new TF1("f1","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *f2  = new TF1("f2","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *fb1 = new TF1("fb1","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *fb2 = new TF1("fb2","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);

  f1->SetParameters(25,5,5);
  f2->SetParameters(25,5,5);
  fb1->SetParameters(25,5,5);
  fb2->SetParameters(25,5,5);
 

  f1->SetParLimits(0,5,50);
  f1->SetParLimits(1,0,20);
  f1->SetParLimits(2,0,20);

  f2->SetParLimits(0,5,150);
  f2->SetParLimits(1,0,20);
  f2->SetParLimits(2,0,20);
  
  fb1->SetParLimits(0,5,50);
  fb1->SetParLimits(1,0,20);
  fb1->SetParLimits(2,0,20);

  fb2->SetParLimits(0,5,150);
  fb2->SetParLimits(1,0,20);
  fb2->SetParLimits(2,0,20);
*/

  pw->GetXaxis()->SetRangeUser(10770,16000);
  pw->GetXaxis()->SetLabelSize(0.08);
  pw->GetYaxis()->SetLabelSize(0.08);
  pw->SetStats(0);



  current = 10770;

  TCanvas *c = new TCanvas;
  c->SetTitle(Form("Rebin = %i",rebin));
  c->Divide(5,4);
  //TCanvas *c1 = new TCanvas;
  //c1->Divide(5,4);

  TGraphErrors *gr = new TGraphErrors(20);
  for(int x=0;x<20;x++) {
    
    //TLine *l5 = new TLine(146,0,146,2000);
    //l5->SetLineColor(kRed); 
    //TLine *l6 = new TLine(152,0,152,2000);
    //l6->SetLineColor(kRed); 
    //TLine *l3 = new TLine(154,0,154,2000);
    //l3->SetLineColor(kBlue); 
    //TLine *l4 = new TLine(160,0,160,2000);
    //l4->SetLineColor(kBlue); 
    //singles[x]->GetXaxis()->SetRangeUser(140,194);
    //c1->cd(x+1);
    //singles[x]->Draw();
    //l3->Draw("same");
    //l4->Draw("same");
    //l5->Draw("same");
    //l6->Draw("same");

    f[x]  = new TF1(Form("f%i",x),"[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
    f[x]->SetParameters(25,5,5);
    f[x]->SetParLimits(0,5,2000);
    f[x]->SetParLimits(1,0,20);
    f[x]->SetParLimits(2,0,20);



    c->cd(x+1);
    pb[x]->Draw();
    pb[x]->SetStats(0);
    pb[x]->Fit(f[x]);

     c->Update();
     double y = gPad->GetUymax();

    TLatex l;
    l.SetTextSize(0.09);
    l.DrawLatex(15,y*.4,Form("t_{1/2} = %f",f[x]->GetParameter(1)));

    TPad *p = new TPad(Form("p%i",x),"",.6,.6,1,1);
    p->Draw();
    p->cd();
    pw->Draw("colz same");   

    c->Update();
 
    TLine *l1=new TLine(current,p->GetUymin(),current,p->GetUymax());
    TLine *l2=new TLine(current+width,p->GetUymin(),current+width,p->GetUymax());
    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l1->Draw();
    l2->Draw();
    current+=shift;
    
    //gr->SetPoint(x,x+1,f[x]->GetParameter(1));
    gr->SetPoint(x,current+width/2,f[x]->GetParameter(1));
    gr->SetPointError(x,0,f[x]->GetParError(1));
 
  }
  new TCanvas;
  gr->Draw("A*");


}
